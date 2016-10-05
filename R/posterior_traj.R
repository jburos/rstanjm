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

#' Estimate the marginal or subject-specific longitudinal trajectory 
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
#' @seealso \code{\link{pp_check}} for graphical posterior predictive checks
#'   of the longitudinal submodel(s).
#'   
#' @examples
#' 
posterior_traj <- function(object, m = 1, newdata = NULL, ids, 
                           limits = c(.025, .975),
                           last_time = NULL, 
                              interpolate = TRUE, extrapolate = TRUE,
                              interpolate_args = list(increments = 25),
                              extrapolate_args = list(dist = NULL, prop = 0.5, 
                                                      increments = 25),
                              re.form = NULL, fun = NULL, 
                              draws = NULL, seed = NULL, offset = NULL, ...) {
  validate_stanjm_object(object)
  M <- object$n_markers
  id_var <- object$id_var
  time_var <- object$time_var
  if (missing(ids)) ids <- NULL
  if (m < 1)
    stop("'m' must be positive")
  if (m > M)
    stop(paste0("'m' must be less than, or equal to, the number ", 
                "of longitudinal markers (M=", M, ")"))
  if (!is.null(seed)) 
    set.seed(seed)
  if (!is.null(fun)) 
    fun <- match.fun(fun)
  
  mf <- model.frame(object)[[m]]
  
  if (!is.null(newdata)) {      
    obsdata <- newdata <- as.data.frame(newdata)
    if (any(is.na(newdata))) 
      stop("Currently NAs are not allowed in 'newdata'.")     
    ids_and_times <- newdata[, c(id_var, time_var)]
    # Latest known observation time for each individual
    if (is.null(last_time))
      last_time <- tapply(ids_and_times[[time_var]], ids_and_times[[id_var]], FUN = max)
  } else {
    # Latest known observation time for each individual
    if (!is.null(last_time))
      stop("'last_time' cannot be specified if 'newdata' is not specified. The ",
           "latest observation (event or censoring) time will be inferred from ",
           "the original data.")
    obsdata <- newdata <- as.data.frame(mf)
    last_time <- object$eventtime
    if (!is.null(ids)) {
      newdata <- newdata[newdata[[id_var]] %in% ids,]
      obsdata <- obsdata[obsdata[[id_var]] %in% ids,]
      last_time <- last_time[as.character(ids)]
    }     
  }
  
  # Time sequence across which to generate the longitudinal trajectory
  id_list <- names(last_time)
  id_class <- class(newdata[[id_var]])
  class(id_list) <- id_class
  max_time <- max(object$eventtime)
  prop <- extrapolate_args$prop
  iinc <- interpolate_args$increments
  einc <- extrapolate_args$increments
  # subject-specific predictions
  if (interpolate) {  
    # Note: prop assumes all individuals entered at time 0
    inter_time_seq <- sapply(0:iinc, function(x, t) t * (x / iinc), 
                       t = last_time)
    # if there is only one patient then need to transform
    if (is.vector(inter_time_seq)) {
      inter_time_seq <- t(inter_time_seq)
      rownames(inter_time_seq) <- id_list
    }
  }  
  if (extrapolate) {  
    # Note: prop assumes all individuals entered at time 0
    dist <- if (!is.null(prop)) prop * (last_time - 0) else
      extrapolate_args$dist
    extra_time_seq <- sapply(0:einc, function(x, t) t + dist * (x / einc), 
                       t = last_time)
    # if there is only one patient then need to transform
    if (is.vector(extra_time_seq)) {
      extra_time_seq <- t(extra_time_seq)
      rownames(extra_time_seq) <- id_list
    }
  }
  if (interpolate && extrapolate) {
    time_seq <- cbind(inter_time_seq, extra_time_seq)
  } else if (interpolate) {
    time_seq <- inter_time_seq
  } else {
    time_seq <- extra_time_seq
  }
  time_seq <- as.data.frame(time_seq)
  n_obs <- NCOL(time_seq)
  colnames(time_seq) <- paste0("V", 1:n_obs)
  time_seq$id <- rownames(time_seq)
  time_seq <- reshape(time_seq, direction = "long", varying = paste0("V", 1:n_obs), 
                 v.names = time_var, timevar = "obs", idvar = id_var)
  class(time_seq[[id_var]]) <- id_class
  class(newdata[[id_var]]) <- id_class
  newdata[[time_var]] <- as.numeric(newdata[[time_var]])
  newdata <- data.table::data.table(newdata, key = c(id_var, time_var))
  newdata <- newdata[data.table::SJ(time_seq[[id_var]], time_seq[[time_var]]),
                         roll = TRUE, rollends = c(TRUE, TRUE)]  
  
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
  ytilde_med <- apply(ytilde, 2, median)
  ytilde_lb <- apply(ytilde, 2, quantile, limits[1]) 
  ytilde_ub <- apply(ytilde, 2, quantile, limits[2]) 

  mutilde <- ppargs$mu
  if (!is.null(newdata) && nrow(newdata) == 1L) 
    mutilde <- t(mutilde)
  mutilde_med <- apply(mutilde, 2, median)
  mutilde_lb <- apply(mutilde, 2, quantile, limits[1]) 
  mutilde_ub <- apply(mutilde, 2, quantile, limits[2])   
    
  Terms <- terms(mf)
  yvar <- rownames(attr(Terms, "factors"))[attr(Terms, "response")]
  xvars <- rownames(attr(Terms, "factors"))[-attr(Terms, "response")]
  dat <- if (length(xvars)) as.data.frame(newdata)[, xvars, drop = FALSE] else NULL
  class(dat[[id_var]]) <- id_class
  out <- data.frame(cbind(dat, ypred = ytilde_med, 
                          ci_lb = mutilde_lb, ci_ub = mutilde_ub,
                          pi_lb = ytilde_lb, pi_ub = ytilde_ub))
  
  class(out) <- c("predict.stanjm", "data.frame")
  structure(out, observed_data = obsdata, 
            y_var = yvar, id_var = id_var, time_var = time_var,
            interpolate = interpolate, interpolate_args = interpolate_args,
            extrapolate = extrapolate, extrapolate_args = extrapolate_args, 
            ids = id_list, last_time = last_time, 
            draws = draws, fun = fun, seed = seed)  
}



