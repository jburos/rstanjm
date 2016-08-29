# Part of the rstanjm package
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

#' Draw from posterior predictive distribution for the longitudinal submodel
#' 
#' The posterior predictive distribution is the distribution of the outcome 
#' implied by the model after using the observed data to update our beliefs 
#' about the unknown parameters in the model. Simulating data from the posterior
#' predictive distribution using the observed predictors is useful for checking 
#' the fit of the model. Drawing from the posterior predictive distribution at 
#' interesting values of the predictors also lets us visualize how a 
#' manipulation of a predictor affects (a function of) the outcome(s). With new 
#' observations of predictor variables we can use the posterior predictive 
#' distribution to generate predicted outcomes.
#' 
#' @export
#' @templateVar stanjmArg object
#' @templateVar mArg object
#' @template args-stanjm-object
#' @template args-m
#' @param newdata Optionally, a data frame in which to look for variables with 
#'   which to predict. If omitted, the model matrix is used. If \code{newdata} 
#'   is provided and any variables were transformed (e.g. rescaled) in the data 
#'   used to fit the model, then these variables must also be transformed in 
#'   \code{newdata}. This only applies if variables were transformed before 
#'   passing the data to one of the modeling functions and \emph{not} if 
#'   transformations were specified inside the model formula. Also see the Note
#'   section below for a note about using the \code{newdata} argument with with
#'   binomial models.
#' @param draws An integer indicating the number of draws to return. The default
#'   and maximum number of draws is the size of the posterior sample.
#' @param re.form If \code{object} contains \code{\link[=stan_glmer]{group-level}}
#'   parameters, a formula indicating which group-level parameters to 
#'   condition on when making predictions. \code{re.form} is specified in the 
#'   same form as for \code{\link[lme4]{predict.merMod}}. The default, 
#'   \code{NULL}, indicates that all estimated group-level parameters are 
#'   conditioned on. To refrain from conditioning on any group-level parameters,
#'   specify \code{NA} or \code{~0}. The \code{newdata} argument may include new
#'   \emph{levels} of the grouping factors that were specified when the model 
#'   was estimated, in which case the resulting posterior predictions 
#'   marginalize over the relevant variables.
#' @param fun An optional function to apply to the results. \code{fun} is found 
#'   by a call to \code{\link{match.fun}} and so can be specified as a function
#'   object, a string naming a function, etc.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param offset A vector of offsets. Only required if \code{newdata} is
#'   specified and an \code{offset} argument was specified when fitting the
#'   model.
#' @param ... Currently unused.
#' 
#' @return A \code{draws} by \code{nrow(newdata)} matrix of simulations
#'   from the posterior predictive distribution. Each row of the matrix is a
#'   vector of predictions generated using a single draw of the model parameters
#'   from the posterior distribution.
#'   
#' @note For binomial models with a number of trials greater than one (i.e., not
#'   Bernoulli models), if \code{newdata} is specified then it must include all 
#'   variables needed for computing the number of binomial trials to use for the
#'   predictions. For example if the left-hand side of the model formula is 
#'   \code{cbind(successes, failures)} then both \code{successes} and 
#'   \code{failures} must be in \code{newdata}. The particular values of 
#'   \code{successes} and \code{failures} in \code{newdata} do not matter so 
#'   long as their sum is the desired number of trials. If the left-hand side of
#'   the model formula were \code{cbind(successes, trials - successes)} then
#'   both \code{trials} and \code{successes} would need to be in \code{newdata},
#'   probably with \code{successes} set to \code{0} and \code{trials} specifying
#'   the number of trials. See the Examples section below and the 
#'   \emph{How to Use the rstanarm Package} for examples.
#' 
#' @seealso \code{\link{pp_check}} for graphical posterior predictive checks.
#'   Examples of posterior predictive checking can also be found in the
#'   \pkg{rstanarm} vignettes and demos.
#'   
#' @examples
#' 
posterior_predict <- function(object, m = 1, newdata = NULL, draws = NULL, 
                              re.form = NULL, fun = NULL, seed = NULL, 
                              offset = NULL, ...) {
  validate_stanjm_object(object)
  M <- object$n_markers
  if (m < 1)
    stop("'m' must be positive")
  if (m > M)
    stop(paste0("'m' must be less than, or equal to, the number ", 
                "of longitudinal markers (M=", M, ")"))
  if (!is.null(seed)) 
    set.seed(seed)
  if (!is.null(fun)) 
    fun <- match.fun(fun)
  if (!is.null(newdata)) {      
    newdata <- as.data.frame(newdata)
    if (any(is.na(newdata))) 
      stop("Currently NAs are not allowed in 'newdata'.")     
  }
            
  dat <-
    pp_data(object,
            m = m,
            newdata = newdata,
            re.form = re.form,
            offset = offset,
            ...)
  ppargs <- pp_args(object, data = pp_eta(object, dat, m, draws), m)
  if (rstanarm:::is.binomial(family(object)[[m]]$family))
    ppargs$trials <- pp_binomial_trials(object, newdata)
  
  ppfun <- pp_fun(object, m)
  ytilde <- do.call(ppfun, ppargs)
  if (!is.null(newdata) && nrow(newdata) == 1L) 
    ytilde <- t(ytilde)
  if (!is.null(fun)) 
    ytilde <- do.call(fun, list(ytilde))
  
  ytilde
}


# call to rstanarm, to get the function to draw from the 
# appropriate posterior predictive distribution for marker m
#
# @param object A stanjm object
# @param m Integer specifying the longitudinal submodel
pp_fun <- function(object, m) {
  suffix <- family(object)[[m]]$family
  getFromNamespace(paste0(".pp_", suffix), "rstanarm")
}


# create list of arguments to pass to the function returned by pp_fun
#
# @param object A stanjm object
# @param data Output from pp_eta (named list with eta and stanmat)
# @return named list
pp_args <- function(object, data, m) {
  stanmat <- data$stanmat
  eta <- data$eta
  stopifnot(is.stanjm(object), is.matrix(stanmat))
  inverse_link <- linkinv(object)[[m]]
  
  args <- list(mu = inverse_link(eta))
  famname <- family(object)[[m]]$family
  if (rstanarm:::is.gaussian(famname)) {
    args$sigma <- stanmat[, paste0("Long", m, "|sigma")]
  } else if (rstanarm:::is.gamma(famname)) {
    args$shape <- stanmat[, paste0("Long", m, "|shape")]
  } else if (rstanarm:::is.ig(famname)) {
    args$lambda <- stanmat[, paste0("Long", m, "|lambda")]
  } else if (rstanarm:::is.nb(famname)) {
    args$size <- stanmat[, paste0("Long", m, "|overdispersion")]
  }
  args
}

# create eta and stanmat (matrix of posterior draws)
# 
# @param object stanjm object
# @param data output from pp_data()
# @param draws number of draws
# @return linear predictor "eta" and matrix of posterior draws stanmat"
pp_eta <- function(object, data, m, draws = NULL) {
  M <- object$n_markers
  x <- data$x
  S <- rstanarm:::posterior_sample_size(object)
  if (is.null(draws)) 
    draws <- S
  if (draws > S) {
    err <- paste0("'draws' should be <= posterior sample size (", 
                  S, ").")
    stop(err)
  }
  some_draws <- isTRUE(draws < S)
  if (some_draws)
    samp <- sample(S, draws)
  if (is.null(data$Zt)) {
    stanmat <- as.matrix.stanjm(object)
    nms <- collect_nms(colnames(stanmat), M)   
    beta <- stanmat[, nms$y[[m]], drop = FALSE]
    if (some_draws) 
      beta <- beta[samp, , drop = FALSE]
    eta <- rstanarm:::linear_predictor.matrix(beta, x, data$offset)
  } else {
    stanmat <- as.matrix(object$stanfit)
    nms <- collect_nms(colnames(stanmat), M)   
    beta <- stanmat[, nms$y[[m]], drop = FALSE]
    if (some_draws) 
      beta <- beta[samp, , drop = FALSE]
    eta <- rstanarm:::linear_predictor.matrix(beta, x, data$offset)
    b <- stanmat[, nms$y_b[[m]], drop = FALSE]
    if (some_draws) 
      b <- b[samp, , drop = FALSE]
    if (is.null(data$Z_names)) {
      b <- b[, !grepl("_NEW_", colnames(b), fixed = TRUE), drop = FALSE]
    } else {
      b <- rstanarm:::pp_b_ord(b, data$Z_names)
    }
    eta <- eta + as.matrix(b %*% data$Zt)
  }
  rstanarm:::nlist(eta, stanmat)
}

pp_b_ord <- function(b, Z_names) {
  ord <- sapply(Z_names, FUN = function(x) {
    m <- grep(paste0("b[", x, "]"), colnames(b), fixed = TRUE)
    len <- length(m)
    if (len == 1) 
      return(m)
    if (len > 1) 
      stop("multiple matches bug")
    m <- grep(paste0("b[", sub(" (.*):.*$", " \\1:_NEW_\\1", x), "]"), 
              colnames(b), fixed = TRUE)
    if (len == 1)
      return(m)
    if (len > 1)
      stop("multiple matches bug")
    x <- strsplit(x, split = ":", fixed = TRUE)[[1]]
    stem <- strsplit(x[[1]], split = " ", fixed = TRUE)[[1]]
    x <- paste(x[1], x[2], paste0("_NEW_", stem[2]), x[2], sep = ":")
    m <- grep(paste0("b[", x, "]"), colnames(b), fixed = TRUE)
    len <- length(m)
    if (len == 1)
      return(m)
    if (len > 1)
      stop("multiple matches bug")
    x <- paste(paste(stem[1], stem[2]), paste0("_NEW_", stem[2]), sep = ":")
    m <- grep(paste0("b[", x, "]"), colnames(b), fixed = TRUE)
    len <- length(m)
    if (len == 1)
      return(m)
    if (len > 1)
      stop("multiple matches bug")
    stop("no matches bug")    
  })
  b[, ord, drop = FALSE]
}

# Number of trials for binomial models
pp_binomial_trials <- function(object, newdata = NULL, m) {
  y <- if (is.null(newdata))
    get_y(object)[[m]] else eval(formula(object)[[m]][[2L]], newdata)
  if (NCOL(y) == 2L) 
    return(rowSums(y))
  rep(1, NROW(y))
}
