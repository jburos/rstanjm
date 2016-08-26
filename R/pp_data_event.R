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


# @param tcalc A vector with the times at which to evaluate the design matrix. The
#   first Npat rows of the design matrix will correspond to time tcalc, with the
#   remaining rows corresponding to design matrix evaluated at the quadrature points
pp_data_event <-
  function(object, newdataEvent, newdataLong,
           offset = NULL,
           ids, times, id_var, time_var,
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
    if (!identical(Npat_event, length(times)))
      stop("Bug found: the vector in the times argument -- which is a vector indicating ",
           "the times at which  to evaluate the survival probability -- should be ",
           "the same length as the number of individuals in newdataEvent.")
    
    # Quadrature points
    quadnodes <- object$quadnodes  # num. of quadrature nodes
    quadpoints <- get_quadpoints(quadnodes) # standardised quadrature points and weights
    quadpoint <- lapply(quadpoints$points, FUN = function(x) (times/2) * x + (times/2))
    quadpoint <- c(list(times), quadpoint)
    
    # Call to evaluate design matrices at quadrature points
    e_xQ <- .pp_data_event_xQ(object, newdataEvent, quadpoint, 
                              ids, id_var, time_var, ...)
    y_xQ <- lapply(1:M, function(m) 
                   .pp_data_long_xQ(object, newdataLong[[m]], m, quadpoint, 
                                    ids, id_var, time_var, ...))
    y_zQ <- lapply(1:M, function(m) 
                   .pp_data_long_zQ(object, newdataLong[[m]], m, quadpoint, 
                                    ids, id_var, time_var, ...))
    #offset <- .pp_data_offset(object, newdata, offset)
    return(rstanarm:::nlist(e_xQ, y_xQ, y_zQ, offset = offset, quadpoint = quadpoint))
  }


# @param newdata A data frame indicating the predictor values for the event
#   submodel. This should include variables id_var and time_var.
# @param quadpoint A list of numeric vectors. The first element of the list is
#   the times at which to evaluate the survival probability, and the remaining
#   elements contain the times for each of the quadrature points.
# @param id_var Name of the id variable
# @param time_var Name of the time variable
.pp_data_event_xQ <- function(object, newdata, quadpoint,
                              ids, id_var, time_var, ...) {
  # Obtain model formula and original factor levels
  form <- formula(object)$Event
  L <- length(form)
  RHS <- formula(substitute(~R, list(R = form[[L]])))
  Terms <- terms(object)$Event
  mf <- model.frame(object)$Event
  ff <- formula(form)
  vars <- rownames(attr(terms.formula(ff), "factors"))
  mf <- mf[vars]
  isFac <- vapply(mf, is.factor, FUN.VALUE = TRUE)
  isFac[attr(Terms, "response")] <- FALSE
  orig_levs <- if (length(isFac) == 0) 
    NULL else lapply(mf[isFac], levels)
  # Evaluate predictor values at quadrature points
  resp_type <- attr(newdata[[1]], "type")
  if (!is.null(resp_type)) {
    newdata <- cbind(unclass(newdata[[1]]), newdata)
    if (resp_type == "right") {
      newdata <- data.table::data.table(newdata, key = c(id_var, "time"))
    } else if (resp_type == "counting") {
      newdata <- data.table::data.table(newdata, key = c(id_var, "start"))
    } else {
      stop("Bug found: newdataEvent was set to NULL, but the model ",
           "frame collected from the original model doesn't appear to ",
           "contain an appropriate Surv(.) response variable.")
    }
  } else {
    newdata <- data.table::data.table(newdata, key = c(id_var, time_var))
  }
  newdataQ <- do.call(rbind, lapply(quadpoint, FUN = function(x) 
    newdata[data.table::SJ(ids, x), roll = TRUE, rollends = c(TRUE, TRUE)])) 
  # Generate model frame based on data at all quadrature points
  mfnew <- model.frame(delete.response(Terms), newdataQ, xlev = orig_levs)
  x <- model.matrix(RHS, data = mfnew, contrasts.arg = attr(mf, "contrasts"))
  return(x)
}

# @param newdata The list or data frame supplying the new longitudinal data to 
#   be used for estimating eta for longitudinal submodel m
# @param m Integer specifying the longitudinal submodel
# @param quadpoint A list of numeric vectors. The first element of the list is
#   the times at which to evaluate the survival probability, and the remaining
#   elements contain the times for each of the quadrature points.
# @param id_var Name of the id variable
# @param time_var Name of the time variable
.pp_data_long_xQ <- function(object, newdata, m, quadpoint, 
                             ids, id_var, time_var, ...) {
  form <- formula(object)[[m]]
  L <- length(form)
  form[[L]] <- lme4::nobars(form[[L]])
  RHS <- formula(substitute(~R, list(R = form[[L]])))
  Terms <- terms(object)[[m]] 
  mf <- model.frame(object)[[m]]
  ff <- formula(form)
  vars <- rownames(attr(terms.formula(ff), "factors"))
  mf <- mf[vars]
  isFac <- vapply(mf, is.factor, FUN.VALUE = TRUE)
  isFac[attr(Terms, "response")] <- FALSE
  orig_levs <- if (length(isFac) == 0) 
    NULL else lapply(mf[isFac], levels)
  # Evaluate predictor values at quadrature points
  newdata <- data.table::data.table(newdata, key = c(id_var, time_var))
  newdataQ <- do.call(rbind, lapply(quadpoint, FUN = function(x) 
    newdata[data.table::SJ(ids, x), roll = TRUE, rollends = c(TRUE, TRUE)]))  
  # Generate model frame based on data at all quadrature points
  mfnew <- model.frame(delete.response(Terms), newdataQ, xlev = orig_levs)
  x <- model.matrix(RHS, data = mfnew, contrasts.arg = attr(mf, "contrasts"))
  return(x)
}

.pp_data_long_zQ <- function(object, newdata, m, quadpoint, 
                             ids, id_var, time_var, re.form = NULL,
                           allow.new.levels = TRUE, na.action = na.pass) {
  # Evaluate predictor values at quadrature points
  newdata <- data.table::data.table(newdata, key = c(id_var, time_var))
  newdata <- do.call(rbind, lapply(quadpoint, FUN = function(x) 
    newdata[data.table::SJ(ids, x), roll = TRUE, rollends = c(TRUE, TRUE)]))  
  # Carry out all remaining steps on expanded newdata
  mfnew <- model.frame(delete.response(terms(object, fixed.only = TRUE)[[m]]),
                       newdata, na.action = na.action)  
  newdata.NA <- newdata
  if (!is.null(fixed.na.action <- attr(mfnew,"na.action"))) {
    newdata.NA <- newdata.NA[-fixed.na.action,]
  }
  tt <- delete.response(terms(object, random.only = TRUE)[[m]])
  rfd <- model.frame(tt, newdata.NA, na.action = na.pass)
  if (!is.null(fixed.na.action))
    attr(rfd,"na.action") <- fixed.na.action
  if (is.null(re.form)) 
    re.form <- rstanarm:::justRE(formula(object)[[m]])
  if (!inherits(re.form, "formula"))
    stop("'re.form' must be NULL, NA, or a formula.")
  if (length(fit.na.action <- attr(mfnew,"na.action")) > 0) {
    newdata <- newdata[-fit.na.action,]
  }
  ReTrms <- lme4::mkReTrms(lme4::findbars(re.form[[2]]), rfd)
  if (!allow.new.levels && any(vapply(ReTrms$flist, anyNA, NA)))
    stop("NAs are not allowed in prediction data",
         " for grouping variables unless 'allow.new.levels' is TRUE.")
  ns.re <- names(re <- ranef(object)[[m]])
  nRnms <- names(Rcnms <- ReTrms$cnms)
  if (!all(nRnms %in% ns.re))
    stop("Grouping factors specified in re.form that were not present in original model.")
  new_levels <- lapply(ReTrms$flist, function(x) levels(factor(x)))
  Zt <- ReTrms$Zt
  p <- sapply(ReTrms$cnms, FUN = length)
  l <- sapply(attr(ReTrms$flist, "assign"), function(i) 
    nlevels(ReTrms$flist[[i]]))
  t <- length(p)
  group_nms <- names(ReTrms$cnms)
  Z_names <- character()
  for (i in seq_along(ReTrms$cnms)) {
    # if you change this, change it in stan_glm.fit() as well
    nm <- group_nms[i]
    nms_i <- paste(ReTrms$cnms[[i]], group_nms[i])
    if (length(nms_i) == 1) {
      Z_names <- c(Z_names, paste0("Long", m, "|", nms_i, ":", levels(ReTrms$flist[[nm]])))
    } else {
      Z_names <- c(Z_names, c(t(sapply("Long", m, "|", nms_i, paste0, ":", new_levels[[nm]]))))
    }
  }
  z <- rstanarm:::nlist(Zt = ReTrms$Zt, Z_names)
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
