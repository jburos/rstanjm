# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
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

#' Methods for stanjm objects
#' 
#' S3 methods for \link[=stanjm-object]{stanjm} objects. There are also 
#' several methods (listed in See Also, below) with their own individual help
#' pages.
#' 
#' @name stanjm-methods
#' @aliases VarCorr fixef ranef ngrps sigma
#' 
#' @templateVar stanjmArg object,x
#' @template args-stanjm-object
#' @param ... Ignored, except by the \code{update} method. See
#'   \code{\link{update}}.
#' 
#' @details Most of these methods are similar to the methods defined for objects
#'   of class 'lm', 'glm', 'glmer', etc. However there are a few exceptions:
#'   
#' \describe{
#' #' \item{\code{coef}}{
#' Returns a list equal to the length of the number of submodels. The first
#' elements of the list are the coefficients from each of the fitted longitudinal
#' submodels and are the same layout as those returned by \code{coef} method of the
#' \pkg{lme4} package, that is, the sum of the random and fixed effects coefficients 
#' for each explanatory variable for each level of each grouping factor. The final
#' element of the returned list is a vector of fixed effect coefficients from the
#' event submodel. 
#' }
#' 
#' \item{\code{confint}}{
#' For models fit using optimization, confidence intervals are returned via a
#' call to \code{\link[stats]{confint.default}}. If \code{algorithm} is
#' \code{"sampling"}, \code{"meanfield"}, or \code{"fullrank"}, the
#' \code{\link{posterior_interval}} function should be used to compute Bayesian
#' uncertainty intervals.
#' }
#' 
#' \item{\code{log_lik}}{
#' For models fit using MCMC only, the \code{log_lik}
#' function returns the \eqn{S} by \eqn{N} pointwise log-likelihood matrix,
#' where \eqn{S} is the size of the posterior sample and \eqn{N} is the number
#' of data points. Note: we use \code{log_lik} rather than defining a
#' \code{\link[stats]{logLik}} method because (in addition to the conceptual
#' difference) the documentation for \code{logLik} states that the return value
#' will be a single number, whereas we return a matrix.
#' }
#' 
#' \item{\code{residuals}}{
#' Residuals are \emph{always} of type \code{"response"} (not \code{"deviance"}
#' residuals or any other type). However, in the case of \code{\link{stan_polr}}
#' with more than two response categories, the residuals are the difference 
#' between the latent utility and its linear predictor.
#' }
#' \item{\code{coef}}{
#' Medians are used for point estimates. See the \emph{Point estimates} section
#' in \code{\link{print.stanjm}} for more details.
#' }
#' \item{\code{se}}{
#' The \code{se} function returns standard errors based on 
#' \code{\link{mad}}. See the \emph{Uncertainty estimates} section in
#' \code{\link{print.stanjm}} for more details.
#' }
#' }
#' 
#' @seealso
#' Other S3 methods for stanjm objects, which have separate documentation, 
#' including \code{\link{as.matrix.stanjm}}, \code{\link{plot.stanjm}}, 
#' \code{\link{predict.stanjm}}, \code{\link{print.stanjm}}, and
#' \code{\link{summary.stanjm}}.
#' 
#' \code{\link{posterior_interval}} and \code{\link{posterior_predict}} for 
#' alternatives to \code{confint} and \code{predict} for models fit using MCMC 
#' or variational approximation.
#' 
NULL


#' @rdname stanjm-methods
#' @export
#' @export coef
#' @param remove_stub Logical specifying whether or not to remove the string
#'    identifying the submodel from each of the coefficient names
#'    
coef.stanjm <- function(object, ...) {
  M <- object$n_markers
  if (length(list(...))) 
    warning("Arguments named \"", paste(names(list(...)), collapse = ", "), 
            "\" ignored.", call. = FALSE)
  fef <- lapply(fixef(object), function(x) data.frame(rbind(x), check.names = FALSE))
  ref <- ranef(object)
  refnames <- lapply(ref, function(x) unlist(lapply(x, colnames)))
  missnames <- lapply(1:M, function(m) setdiff(refnames[[m]], names(fef[[m]])))
  nmiss <- sapply(missnames, length)
  if (any(nmiss > 0)) for (m in 1:M) {
    if (nmiss[m] > 0) {
      fillvars <- setNames(data.frame(rbind(rep(0, nmiss[m]))), missnames[[m]])
      fef[[m]] <- cbind(fillvars, fef[[m]])
    }
  }
  val <- lapply(1:M, function(m) 
    lapply(ref[[m]], function(x) fef[[m]][rep.int(1L, nrow(x)), , drop = FALSE]))
  for (m in 1:M) {  # loop over number of markers
    for (i in seq(a = val[[m]])) {  # loop over number of grouping factors
      refi <- ref[[m]][[i]]
      row.names(val[[m]][[i]]) <- row.names(refi)
      nmsi <- colnames(refi)
      if (!all(nmsi %in% names(fef[[m]]))) 
        stop("Unable to align random and fixed effects.", call. = FALSE)
      for (nm in nmsi) 
        val[[m]][[i]][[nm]] <- val[[m]][[i]][[nm]] + refi[, nm]
    }
  }
  val <- lapply(val, function(x) structure(x, class = "coef.mer"))
  val <- c(val, fef[length(fef)])         
  val <- list_nms(val, M)
  val        
}

#' @rdname stanjm-methods
#' @export
#' 
fitted.stanjm <- function(object, ...)  {
  object$fitted.values
}

#' Extract pointwise log-likelihood matrix
#' 
#' @export
#' @keywords internal
#' @param object Fitted model object.
#' @param ... Arguments to methods. For example the
#'   \code{\link[=stanjm-methods]{stanjm}} method accepts the argument
#'   \code{newdata}.
#' @return Pointwise log-likelihood matrix.
#' @seealso \code{\link{log_lik.stanjm}}
#' 
log_lik <- function(object, ...) UseMethod("log_lik")

#' @rdname stanjm-methods
#' @export
#' @param newdata For \code{log_lik}, an optional data frame of new data (e.g. 
#'   holdout data) to use when evaluating the log-likelihood. See the 
#'   description of \code{newdata} for \code{\link{posterior_predict}}.
log_lik.stanjm <- function(object, newdata = NULL, ...) {
  if (!rstanarm:::used.sampling(object)) 
    rstanarm:::STOP_sampling_only("Pointwise log-likelihood matrix")
  if (!is.null(newdata)) {
    newdata <- as.data.frame(newdata)
  }
  fun <- ll_fun(object)
  args <- ll_args(object, newdata)
  sapply(seq_len(args$N), function(i) {
    as.vector(fun(i = i, data = args$data[i, , drop = FALSE], 
                  draws = args$draws))
  })
}

#' @rdname stanjm-methods
#' @export 
residuals.stanjm <- function(object, ...) {
  object$residuals
}

#' Extract standard errors
#' 
#' Generic function for extracting standard errors from a fitted joint model.
#' 
#' @export
#' @keywords internal
#' @param object A fitted model object.
#' @param ... Arguments to methods.
#' @return Standard errors of model parameters.
#' @seealso \code{\link{se.stanjm}}
#' 
se <- function(object, ...) UseMethod("se")

#' @rdname stanjm-methods
#' @export
se.stanjm <- function(object, ...) {
  object$ses
}

#' @rdname stanjm-methods
#' @export
#' @param fixed.only A logical specifying whether to only retain the fixed effect
#'   part of the longitudinal submodel formulas
#' @param random.only A logical specifying whether to only retain the random effect
#'   part of the longitudinal submodel formulas  
formula.stanjm <- function (x, fixed.only = FALSE, random.only = FALSE, ...) {
  if (missing(fixed.only) && random.only) 
    fixed.only <- FALSE
  if (fixed.only && random.only) 
    stop("'fixed.only' and 'random.only' can't both be TRUE.", call. = FALSE)
  
  M <- x$n_markers
  fr <- lapply(x$glmod, function(m) m$fr) 
  form <- lapply(x$formula, as.formula, ...)
  glmod_form <- lapply(seq(M), function(m) attr(fr[[m]], "formula")) 
  if (!is.null(glmod_form)) form[1:M] <- glmod_form[1:M]
  if (any(fixed.only, random.only)) {
    if (fixed.only) {
      for (m in 1:M)
        form[[m]][[length(form[[m]])]] <- lme4::nobars(form[[m]][[length(form[[m]])]])
    }
    if (random.only)
      for (m in 1:M)
        form[[m]] <- rstanarm:::justRE(form[[m]], response = TRUE)
  } 
  return(form)
}

#' terms method for stanjm objects
#' @export
#' @keywords internal
#' @param x,fixed.only,random.only,... See lme4:::terms.merMod.
#' 
terms.stanjm <- function(x, ..., fixed.only = TRUE, random.only = FALSE) {
  if (!is.stanjm(x))
    return(NextMethod("terms"))
  
  fr <- lapply(x$glmod, function(m) m$fr) 
  if (missing(fixed.only) && random.only) 
    fixed.only <- FALSE
  if (fixed.only && random.only) 
    stop("'fixed.only' and 'random.only' can't both be TRUE.", call. = FALSE)
  
  Terms <- lapply(fr, function(m) attr(m, "terms"))
  if (fixed.only) {
    Terms <- lapply(seq(M), function(m) {
      Terms <- terms.formula(formula(x, fixed.only = TRUE)[[m]])
      attr(Terms, "predvars") <- attr(terms(fr[[m]]), "predvars.fixed")
      Terms
    })     
  } 
  if (random.only) {
    Terms <- lapply(seq(M), function(m) {
      Terms <- terms.formula(lme4::subbars(formula.stanjm(x, random.only = TRUE)[[m]]))
      attr(Terms, "predvars") <- attr(terms(fr[[m]]), "predvars.random")
      Terms
    })      
  }
  
  return(Terms)
}


#' @rdname stanjm-methods
#' @export
#' @method update stanjm
#' @param formula.,evaluate See \code{\link[stats]{update}}.
#'
update.stanjm <- function(object, formulaLong., formulaEvent., ..., evaluate = TRUE) {
  call <- getCall(object)
  M <- object$n_markers
  if (is.null(call)) 
    stop("'object' does not contain a 'call' component.", call. = FALSE)
  extras <- match.call(expand.dots = FALSE)$...
  fm <- formula(object)
  if (!missing(formulaLong.)) {
    if (M > 1) {
      if (!is.list(formulaLong.))
        stop("To update the formula for a multivariate joint model ",
             "'formulaLong.' should be a list of formula objects. Use ",
             "'~ .' if you do not wish to alter the formula for one or ",
             "more of the longitudinal submodels.", call. = FALSE)
      if (length(formulaLong.) != M)
        stop(paste0("The list provided in 'formulaLong.' appears to be the ",
             "incorrect length; should be length ", M), call. = FALSE)     
    } else {
      if (!is.list(formulaLong.)) formulaLong. <- list(formulaLong.)
    }
    fm_long <- lapply(1:M, function(m) 
      update.formula(fm[[m]], formulaLong.[[m]]))
    names(fm_long) <- NULL
    fm_long <- as.call(c(quote(list), fm_long))
    call$formulaLong <- fm_long
  }
  if (!missing(formulaEvent.))
    call$formulaEvent <- update.formula(fm[[length(fm)]], formulaEvent.)  
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) 
      call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  
  if (!evaluate) 
    return(call)
  
  # do this like lme4 update.merMod instead of update.default
  ff <- environment(formula(object))
  pf <- parent.frame()
  sf <- sys.frames()[[1L]]
  tryCatch(eval(call, envir = ff),
           error = function(e) {
             tryCatch(eval(call, envir = sf),
                      error = function(e) {
                        eval(call, pf)
                      })
           })
}

#' @rdname stanjm-methods
#' @export 
#' @param correlation For \code{vcov}, if \code{FALSE} (the default) the
#'   covariance matrix is returned. If \code{TRUE}, the correlation matrix is
#'   returned instead.
#'
vcov.stanjm <- function(object, correlation = FALSE, ...) {
  out <- object$covmat
  if (!correlation) return(out)
  out <- lapply(out, cov2cor)
  out
}


.stanjm_check <- function(object) {
  if (!is.stanjm(object))
    stop("This method is for stanjm objects only.", call. = FALSE)
}
.cnms <- function(object, remove_stub = FALSE) {
  .stanjm_check(object)
  cnms <- object$cnms
  if (remove_stub) {
    cnms <- lapply(cnms, rm_stub)
  }
  cnms
}
.p <- function(object) {
  .stanjm_check(object)
  sapply(object$cnms, length)
}
#.y_cnms <- function(object) {
#  .stanjm_check(object)
#  object$y_cnms
#}
#.y_flist <- function(object) {
#  .stanjm_check(object)
#  as.list(object$y_flist)
#}


#' @rdname stanjm-methods
#' @export
#' @export fixef
#' @importFrom lme4 fixef
#' @param remove_stub Logical specifying whether to remove the string identifying 
#'    the longitudinal or event submodel from each of the coefficient names
#' 
fixef.stanjm <- function(object, remove_stub = TRUE, ...) {
  coefs <- object$coefficients
  coefs <- lapply(coefs, function(x) x[b_names(names(x), invert = TRUE)])
  if (remove_stub) {
    for (i in 1:length(coefs)) names(coefs[[i]]) <- rm_stub(names(coefs[[i]]))
  }
  coefs
}

#' @rdname stanjm-methods
#' @export
#' @export ngrps
#' @importFrom lme4 ngrps
#' 
ngrps.stanjm <- function(object, ...) {
  object$n_grps  
}

#' @rdname stanjm-methods
#' @export
#' @export ranef
#' @importFrom lme4 ranef
#'
ranef.stanjm <- function(object, ...) {
  M <- object$n_markers
  all_names <- if (rstanarm:::used.optimizing(object))
    rownames(object$stan_summary) else object$stanfit@sim$fnames_oi
  ans_list <- lapply(1:M, function(m) { 
    sel <- b_names(all_names, m)
    ans <- object$stan_summary[sel, rstanarm:::select_median(object$algorithm)]
    # avoid returning the extra levels that were included
    ans <- ans[!grepl("_NEW_", names(ans), fixed = TRUE)]
    fl <- as.list(object$y_flist[[m]])
    levs <- lapply(fl, levels)
    asgn <- attr(fl, "assign")
    cnms <- object$y_cnms[[m]]
    nc <- vapply(cnms, length, 1L)
    nb <- nc * vapply(levs, length, 1L)[asgn]
    nbseq <- rep.int(seq_along(nb), nb)
    ml <- split(ans, nbseq)
    for (i in seq_along(ml)) {
      ml[[i]] <- matrix(ml[[i]], ncol = nc[i], byrow = TRUE, 
                        dimnames = list(NULL, cnms[[i]]))
    }
    ans <- lapply(seq_along(fl), function(i) {
      data.frame(do.call(cbind, ml[asgn == i]), row.names = levs[[i]], 
                 check.names = FALSE)
    })
    names(ans) <- names(fl)
    structure(ans, class = "ranef.mer")
    ans
  })
  ans_list <- list_nms(ans_list, M)
  ans_list
}


#' @rdname stanjm-methods
#' @export
#' @export sigma
#' @rawNamespace if(getRversion()>='3.3.0') importFrom(stats, sigma) else
#'   importFrom(lme4,sigma)
#'
sigma.stanjm <- function(object, ...) {
  M <- object$n_markers
  nms <- grep("^Long[1-9]\\|sigma", rownames(object$stan_summary), value = TRUE)
  if (!length(nms)) 
    return(1)
  sigma <- object$stan_summary[nms, rstanarm:::select_median(object$algorithm)]
  if (M > 1L) {
    new_nms <- gsub("\\|sigma", "", nms)
    names(sigma) <- new_nms
  }
  return(sigma)
}

#' @rdname stanjm-methods
#' @param sigma Ignored (included for compatibility with
#'   \code{\link[nlme]{VarCorr}}).
#' @export
#' @export VarCorr
#' @importFrom nlme VarCorr
#' @importFrom lme4 mkVarCorr
VarCorr.stanjm <- function(x, sigma = 1, ...) {
  cnms <- .cnms(x)
  means <- rstan::get_posterior_mean(x$stanfit)
  means <- means[, ncol(means)]
  theta <- means[grepl("^theta_L", names(means))]
  p <- .p(x)
  out <- lme4::mkVarCorr(sc = 1, cnms = cnms, 
                         nc = p,
                         theta = theta, nms = names(cnms))
  structure(out, useSc = 0, class = "VarCorr.merMod")
}


# Exported but doc kept internal ----------------------------------------------

#' family method for stanjm objects
#'
#' @keywords internal
#' @export
#' @param object,... See \code{\link[stats]{family}}.
family.stanjm <- function(object, ...) object$family

#' model.frame method for stanjm objects
#' 
#' @keywords internal
#' @export
#' @param formula,... See \code{\link[stats]{model.frame}}.
#' @param fixed.only See \code{\link[lme4]{model.frame.merMod}}.
#' 
model.frame.stanjm <- function(formula, fixed.only = FALSE, ...) {
  if (is.stanjm(formula)) {
    M <- formula$n_markers
    fr <- lapply(seq(M), function(m) formula$glmod[[m]]$fr)
    if (fixed.only) {
      fr <- lapply(seq(M), function(m) {
        ff <- formula(formula, fixed.only = TRUE)[[m]]
        vars <- rownames(attr(terms.formula(ff), "factors"))
        fr[[m]][vars]
      })
    }
    return(fr)
  } 
  
  NextMethod("model.frame")
}


#----------------------------------------------
# The following are not yet adapted for stan_jm
#----------------------------------------------

#' model.matrix method for stanjm objects
#' 
#' @keywords internal
#' @export
#' @param object,... See \code{\link[stats]{model.matrix}}.
#' 
model.matrix.stanjm <- function(object, ...) {
  if (is.stanjm(object)) {
    M <- object$n_markers
    return(lapply(seq(M), object$glmod[[m]]$X))
  }
  
  NextMethod("model.matrix")
}

#' terms method for stanjm objects
#' @export
#' @keywords internal
#' @param x,fixed.only,random.only,... See lme4:::terms.merMod.
#' 
terms.stanjm <- function(x, ..., fixed.only = TRUE, random.only = FALSE) {
  if (!is.mer(x))
    return(NextMethod("terms"))
  
  fr <- x$glmod$fr
  if (missing(fixed.only) && random.only) 
    fixed.only <- FALSE
  if (fixed.only && random.only) 
    stop("'fixed.only' and 'random.only' can't both be TRUE.", call. = FALSE)
  
  Terms <- attr(fr, "terms")
  if (fixed.only) {
    Terms <- terms.formula(formula(x, fixed.only = TRUE))
    attr(Terms, "predvars") <- attr(terms(fr), "predvars.fixed")
  } 
  if (random.only) {
    Terms <- terms.formula(lme4::subbars(formula.stanjm(x, random.only = TRUE)))
    attr(Terms, "predvars") <- attr(terms(fr), "predvars.random")
  }
  
  return(Terms)
}
