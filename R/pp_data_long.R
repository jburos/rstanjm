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

pp_data_long <-
  function(object, 
           m,
           newdata = NULL,
           re.form = NULL,
           offset = NULL,
           ...) {
    validate_stanjm_object(object)
    x <- .pp_data_long_x(object, newdata, m, ...)
    z <- .pp_data_long_z(object, newdata, m, re.form, ...)
    #offset <- .pp_data_offset(object, newdata, offset)
    return(rstanarm:::nlist(x, offset = offset, Zt = z$Zt, Z_names = z$Z_names))
  }


# the functions below are heavily based on a combination of 
# lme4:::predict.merMod and lme4:::mkNewReTrms, although they do also have 
# substantial modifications
.pp_data_long_x <- function(object, newdata, m, ...) {
  x <- get_x(object)[[m]]
  if (is.null(newdata)) return(x)
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
  mfnew <- model.frame(delete.response(Terms), newdata, xlev = orig_levs)
  x <- model.matrix(RHS, data = mfnew, contrasts.arg = attr(x, "contrasts"))
  return(x)
}

.pp_data_long_z <- function(object, newdata, m, re.form = NULL,
                           allow.new.levels = TRUE, na.action = na.pass) {
  NAcheck <- !is.null(re.form) && !is(re.form, "formula") && is.na(re.form)
  fmla0check <- (is(re.form, "formula") && 
                   length(re.form) == 2 && 
                   identical(re.form[[2]], 0))
  if (NAcheck || fmla0check) return(list())
  if (is.null(newdata) && is.null(re.form)) {
    Z <- get_z(object)[[m]]
    return(list(Zt = t(Z)))
  }
  else if (is.null(newdata)) {
    rfd <- mfnew <- model.frame(object)[[m]]
  } else {
    if ("gam" %in% names(object))
      stop("'posterior_predict' with non-NULL 're.form' not yet supported ", 
           "for models estimated via 'stan_gamm4'")
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
  }
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
