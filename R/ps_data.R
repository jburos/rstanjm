# Part of the rstanjm package
# Copyright 2015 Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker
# Copyright (C) 2015, 2016 Trustees of Columbia University
# Copyright (C) 2016 Sam Brilleman
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


# @param object A fitted stanjm model
# @param ids A vector of ids indicating the values of the id_var to provide the 
#   predictions for.
# @param t A vector with the times at which to evaluate the survival probability.
# @param id_var Name of the id variable in the original model.
# @param time_var Name of the time variable in the original model.
ps_data <-
  function(object, newdataEvent, newdataLong,
           offset = NULL,
           ids, t, id_var, time_var,
           ...) {
    validate_stanjm_object(object)
    M <- object$n_markers
    flist_long <- lapply(newdataLong, function(x) x[[id_var]])
    flist_event <- newdataEvent[[id_var]]
    Npat_long <- lapply(flist_long, function(x) length(unique(x)))
    Npat_event <- length(unique(flist_event))
    lapply(Npat_long, function(x) {
      if (!identical(x, Npat_event))
        stop("Bug found: number of individuals in newdataEvent and each newdataLong ",
             "should be equal.")})
    if (!identical(Npat_event, length(ids)))
      stop("Bug found: number of individuals in newdataEvent and each newdataLong ",
           "should be the same as the length of the 'ids' argument.")
    if (!identical(Npat_event, length(t)))
      stop("Bug found: the vector in the argument 't' -- which is a vector indicating ",
           "the times at which  to evaluate the survival probability -- should be ",
           "the same length as the number of individuals in newdataEvent.")
    
    # Quadrature points
    quadnodes <- object$quadnodes  # num. of quadrature nodes
    quadpoints <- get_quadpoints(quadnodes) # standardised quadrature points and weights
    tQ <- lapply(quadpoints$points, FUN = function(x) (t/2) * x + (t/2))
    t_and_tQ <- c(list(t), tQ)
    
    # Call to evaluate design matrices at quadrature points
    e_xQ <- .ps_data_event_xQ(object, newdataEvent, t_and_tQ, 
                              ids, id_var, time_var, ...)
    y_xQ <- lapply(1:M, function(m) 
                   .ps_data_long_xQ(object, newdataLong[[m]], m, t_and_tQ, 
                                    ids, id_var, time_var, ...))
    y_zQ <- lapply(1:M, function(m) 
                   .ps_data_long_zQ(object, newdataLong[[m]], m, t_and_tQ, 
                                    ids, id_var, time_var, ...))
    #offset <- .pp_data_offset(object, newdata, offset)
    return(nlist(e_xQ, y_xQ, y_zQ, offset, t, tQ))
  }


# @param object A fitted stanjm model
# @param newdata A data frame indicating the predictor values for the event
#   submodel. This should include variables id_var and time_var.
# @param t_and_tQ A list of numeric vectors. The first element of the list is
#   the times at which to evaluate the survival probability, and the remaining
#   elements contain the times at each of the quadrature points.
# @param ids A vector of ids indicating the values of the id_var to provide the 
#   predictions for.
# @param id_var Name of the id variable in the original model.
# @param time_var Name of the time variable in the original model.
.ps_data_event_xQ <- function(object, newdata, t_and_tQ,
                              ids, id_var, time_var, ...) {
  # Obtain model formula RHS
  form <- formula(object)$Event
  L <- length(form)
  RHS <- formula(substitute(~R, list(R = form[[L]])))
  # Terms object to be used for obtaining distinct
  # variables in the model frame from the fitted model
  Terms <- terms(object)$Event
  mf <- model.frame.stanjm(object)$Event
  ff <- formula(formula(object)$Event)
  vars <- rownames(attr(terms.formula(ff), "factors"))
  # Limit model frame to main effects only and identify
  # which variables are factors
  mf <- mf[vars]
  isFac <- vapply(mf, is.factor, FUN.VALUE = TRUE)
  isFac[attr(Terms, "response")] <- FALSE
  # Obtain the original levels for the factors, to be 
  # used in making sure dumming coding of design matrix
  # based on new data is appropriate
  orig_levs <- if (length(isFac) == 0) 
    NULL else lapply(mf[isFac], levels)
  # Assess whether Surv() response was counting or 
  # right censored -- returns NULL if newdata was specified
  # by the user and not obtained from the fitted model
  resp_type <- attr(newdata[[1]], "type")
  if (!is.null(resp_type)) {
    # newdata was obtained from the fitted model, and therefore time
    # variable taken to be either only observation time (single row
    # per individual) or "start" of time (for multiple row per individual)
    newdata <- cbind(unclass(newdata[[1]]), newdata[,-1])
    if (resp_type == "right") {
      newdata <- data.table::data.table(newdata, key = c(id_var, "time"))
      newdata[["time"]] <- as.numeric(newdata[["time"]])
    } else if (resp_type == "counting") {
      newdata <- data.table::data.table(newdata, key = c(id_var, "start"))
      newdata[["start"]] <- as.numeric(newdata[["start"]])
    } else {
      stop("Bug found: newdataEvent was set to NULL, but the model ",
           "frame collected from the original model doesn't appear to ",
           "contain an appropriate time variable in the Surv(.) response.")
    }
  } else {  
    # newdata was specified by user, and therefore must have included 
    # time_var as one of the variables in the data frame
    if (!time_var %in% colnames(newdata)) 
      newdata[[time_var]] <- rep(0, nrow(newdata))
    newdata <- data.table::data.table(newdata, key = c(id_var, time_var))
  }
  # Expand newdata based on a rolling merge between newdata at observation
  # times and the identified quadrature points
  newdataQ <- lapply(t_and_tQ, function(x) 
    newdata[data.table::SJ(ids, x), roll = TRUE, rollends = c(TRUE, TRUE)]) 
  # Generate new model frames based on newdata at each of the quadrature points
  mfnew <- lapply(newdataQ, function(x) 
    model.frame(delete.response(Terms), x, xlev = orig_levs))
  # Calculate design matrix based on new model frames and using contrasts
  # from the original fitted model
  form_pred <- use_predvars(object$coxmod)
  RHS_pred <- formula(substitute(~R, list(R = form_pred[[L]])))
  x <- lapply(mfnew, function(x) 
    model.matrix(RHS, data = x, contrasts.arg = attr(mf, "contrasts")))
  return(x)
}

# @param object A fitted stanjm model
# @param newdata The list or data frame supplying the new longitudinal data to 
#   be used for estimating eta for longitudinal submodel m
# @param m Integer specifying the longitudinal submodel
# @param t_and_tQ A list of numeric vectors. The first element of the list is
#   the times at which to evaluate the survival probability, and the remaining
#   elements contain the times for each of the quadrature points.
# @param id_var Name of the id variable
# @param time_var Name of the time variable
.ps_data_long_xQ <- function(object, newdata, m, t_and_tQ, 
                             ids, id_var, time_var, ...) {
  # Obtain model formula RHS fixed part only
  form <- formula(object)[[m]]
  L <- length(form)
  form[[L]] <- lme4::nobars(form[[L]])
  RHS <- formula(substitute(~R, list(R = form[[L]])))
  # Terms object to be used for obtaining distinct
  # variables in the model frame from the fitted model
  Terms <- terms(object)[[m]]
  mf <- model.frame.stanjm(object)[[m]]
  ff <- formula(form)
  vars <- rownames(attr(terms.formula(ff), "factors"))
  # Limit model frame to main effects only and identify
  # which variables are factors
  mf <- mf[vars]
  isFac <- vapply(mf, is.factor, FUN.VALUE = TRUE)
  isFac[attr(Terms, "response")] <- FALSE
  # Obtain the original levels for the factors, to be 
  # used in making sure dumming coding of design matrix
  # based on new data is appropriate
  orig_levs <- if (length(isFac) == 0) 
    NULL else lapply(mf[isFac], levels)
  # Expand newdata based on a rolling merge between newdata at observation
  # times and the identified quadrature points
  newdata <- data.table::data.table(newdata, key = c(id_var, time_var))
  newdataQ <- lapply(t_and_tQ, function(x) 
    newdata[data.table::SJ(ids, x), roll = TRUE, rollends = c(TRUE, TRUE)]) 
  # Generate new model frames based on newdata at each of the quadrature points
  mfnew <- lapply(newdataQ, function(x) 
    model.frame(delete.response(Terms), x, xlev = orig_levs))
  # Calculate design matrix based on new model frames and using contrasts
  # from the original fitted model
  form_pred <- use_predvars(object$glmod[[m]])
  form_pred[[L]] <- lme4::nobars(form_pred[[L]])
  RHS_pred <- formula(substitute(~R, list(R = form_pred[[L]])))
  x <- lapply(mfnew, function(x) 
    model.matrix(RHS_pred, data = x, contrasts.arg = attr(mf, "contrasts")))
  return(x)
}

.ps_data_long_zQ <- function(object, newdata, m, t_and_tQ, 
                             ids, id_var, time_var, re.form = NULL,
                           allow.new.levels = TRUE, na.action = na.pass) {
  # Expand newdata based on a rolling merge between newdata at observation
  # times and the identified quadrature points
  newdata <- data.table::data.table(newdata, key = c(id_var, time_var))
  newdataQ <- lapply(t_and_tQ, FUN = function(x) 
    newdata[data.table::SJ(ids, x), roll = TRUE, rollends = c(TRUE, TRUE)])  
  # Evaluate model frame at each quadrature point based on fixed part of
  # original model formula. NAs are allowed to enter into the new model frame
  ttf <- delete.response(terms(object, fixed.only = TRUE)[[m]])
  mfnew <- lapply(newdataQ, function(x)
    model.frame(ttf, x, na.action = na.action))  
  # Identify which rows of the new model frame contain NA and drop
  # those rows from newdataQ.NA
  fixed.na.action <- lapply(mfnew, attr, "na.action")
  newdataQ.NA <- mapply(function(x, y) {
    if (!is.null(y)) x[-y,] else x}, newdataQ, fixed.na.action, SIMPLIFY = FALSE)
  # Evaluate model frame at each quadrature point based on random part of
  # original model formula. (NAs in the random effects part of the model
  # are allowed to pass, but NAs in fixed part were determined by NA action
  # for mfnew).
  ttr <- delete.response(terms(object, random.only = TRUE)[[m]])
  rfd <- lapply(newdataQ.NA, function(x)
    model.frame(ttr, x, na.action = na.pass))
  for (i in 1:length(rfd)) {
    if (!is.null(fixed.na.action[[i]]))
      attr(rfd[[i]],"na.action") <- fixed.na.action[[i]]    
  }
  
  # Get model formula for random effects part
  if (is.null(re.form)) 
    re.form <- justRE(use_predvars(object$glmod[[m]]))
  if (!inherits(re.form, "formula"))
    stop("'re.form' must be NULL, NA, or a formula.")
  
  # Make random effects component
  ReTrms <- lapply(rfd, function(x) {
    ReTrms_tmp <- lme4::mkReTrms(lme4::findbars(re.form[[2]]), x)
    if (!allow.new.levels && any(vapply(ReTrms_tmp$flist, anyNA, NA)))
      stop("NAs are not allowed in prediction data",
           " for grouping variables unless 'allow.new.levels' is TRUE.")
    ReTrms_tmp
  })
  # Component names from original model
  ns.re <- names(re <- ranef(object)[[m]])
  # Component names from ReTrms based on newdata
  nRnms <- sapply(ReTrms, function(x) names(Rcnms <- x$cnms))
  if (!all(nRnms %in% ns.re))
    stop("Grouping factors specified in re.form that were not present in original model.")
  
  # New design matrix and coefficient names
  Zt <- lapply(ReTrms, function(x) x$Zt)
  # List with levels (e.g. ids) for each grouping factor in newdata
  # used for generating Z_names below
  new_levels <- lapply(ReTrms[[1]]$flist, function(x) levels(factor(x)))
  p <- sapply(ReTrms[[1]]$cnms, FUN = length)
  l <- sapply(attr(ReTrms[[1]]$flist, "assign"), function(i) 
    nlevels(ReTrms[[1]]$flist[[i]]))
  t <- length(p)
  group_nms <- names(ReTrms[[1]]$cnms)
  Z_names <- character()
  for (i in seq_along(ReTrms[[1]]$cnms)) {
    # if you change this, change it in stan_glm.fit() as well
    nm <- group_nms[i]
    nms_i <- paste(ReTrms[[1]]$cnms[[i]], group_nms[i])
    if (length(nms_i) == 1) {
      Z_names <- c(Z_names, paste0("Long", m, "|", nms_i, ":", new_levels[[nm]]))
    } else {
      Z_names <- c(Z_names, c(t(sapply(nms_i, function(x) 
        paste0("Long", m, "|", x, ":", new_levels[[nm]])))))
    }
  }

  z <- nlist(Zt, Z_names)
  return(z)
}



# handle offsets ----------------------------------------------------------
null_or_zero <- function(x) {
  isTRUE(is.null(x) || all(x == 0))
}

.pp_data_offset <- function(object, newdata = NULL, offset = NULL) {
  if (is.null(newdata)) {
    # get offset from model object (should be null if no offset)
    if (is.null(offset)) 
      offset <- object$offset %ORifNULL% model.offset(model.frame(object))
  } else {
    if (!is.null(offset))
      stopifnot(length(offset) == nrow(newdata))
    else {
      # if newdata specified but not offset then confirm that model wasn't fit
      # with an offset (warning, not error)
      if (!is.null(object$call$offset) || 
          !null_or_zero(object$offset) || 
          !null_or_zero(model.offset(model.frame(object)))) {
        warning(
          "'offset' argument is NULL but it looks like you estimated ", 
          "the model using an offset term.", 
          call. = FALSE
        )
      }
      offset <- rep(0, nrow(newdata))
    }
  }
  return(offset)
}
