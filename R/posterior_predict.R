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
#' about the unknown parameters in the model. This function allows us to 
#' generate an estimated longitudinal trajectory (either subject-specific, or
#' by marginalising over the distribution of the random effects) based on
#' draws from the posterior predictive distribution.
#' 
#' @export
#' @templateVar stanjmArg object
#' @templateVar mArg m
#' @template args-stanjm-object
#' @template args-m
#' 
#' @param newdata Optionally, a data frame in which to look for variables with 
#'   which to predict. If omitted, the model matrix is used. See the 
#'   \strong{Details} section for further important details that need to 
#'   be considered when specifying \code{newdata}.
#' @param interpolate A logical specifying whether to interpolate the estimated 
#'   longitudinal trajectory in between the observation times. This can be used
#'   to achieve a smooth estimate of the longitudinal trajectory across the 
#'   entire follow up time. If \code{TRUE} then the interpolation can be further 
#'   controlled using the \code{control} argument.
#' @param extrapolate A logical specifying whether to extrapolate the estimated 
#'   longitudinal trajectory beyond the time of the last known observation time.
#'   If \code{TRUE} then the extrapolation can be further controlled using
#'   the \code{control} argument.
#' @param control A named list with parameters controlling the interpolation or
#'   extrapolation of the estimated longitudinal trajectory when either 
#'   \code{interpolate = TRUE} or \code{extrapolate = TRUE}. The 
#'   list can contain one or more of the following named elements: \cr
#'   \describe{
#'     \item{\code{int_points}}{a positive integer specifying the number of discrete 
#'     time points at which to calculate the estimated longitudinal response for
#'     \code{interpolate = TRUE}. These time points are evenly spaced starting at 
#'     0 and ending at the last known observation time for each individual. The
#'     last observation time for each individual is taken to be either: the
#'     event or censoring time if no new data is provided; the time specified
#'     in the "last_time" column if provided in the new data (see \strong{Details}
#'     section below); or the time of the last longitudinal measurement if new
#'     data is provided but no "last_time" column is included. The default is 15.}
#'     \item{\code{ext_points}}{a positive integer specifying the number of discrete 
#'     time points at which to calculate the estimated longitudinal response for
#'     \code{extrapolate = TRUE}. These time points are evenly spaced between the 
#'     last known observation time for each individual and the extrapolation 
#'     distance specifed using either \code{ext_distance} or \code{ext_prop}.
#'     The default is 15.}
#'     \item{\code{ext_prop}}{a positive scalar between 0 and 1 specifying the 
#'     amount of time across which to extrapolate the longitudinal trajectory,
#'     represented as a proportion of the total observed follow up time for each
#'     individual. For example specifying \code{ext_prop = 0.2} means that for an
#'     individual for whom the latest of their measurement, event or censoring times
#'     was 10 years, their estimated longitudinal trajectory will be extrapolated 
#'     out to 12 years (i.e. 10 + (0.2 * 10)). The default value is 0.2.}
#'     \item{\code{ext_distance}}{a positive scalar specifying the amount of time 
#'     across which to extrapolate the longitudinal trajectory for each individual,
#'     represented in units of the time variable \code{time_var} (from fitting the
#'     model). This cannot be specified if \code{ext_prop} is specified.} 
#'   }
#' @param ids An optional vector specifying a subset of IDs for whom the 
#'   predictions should be obtained. The default is to predict for all individuals
#'   who were used in estimating the model or, if \code{newdata} is specified,
#'   then all individuals contained in \code{newdata}.
#' @param prob A scalar between 0 and 1 specifying the width to use for the 
#'   uncertainty interval (sometimes called credible interval) for the predicted
#'   mean response and the prediction interval for the predicted (raw) response. 
#'   For example \code{prob = 0.95} (the default) means that the 2.5th and 97.5th  
#'   percentiles will be provided.
#' @param re.form A formula indicating which random effects (group-level parameters)
#'   to condition on when making predictions. \code{re.form} is specified in the 
#'   same form as for \code{\link[lme4]{predict.merMod}}. The default, 
#'   \code{NULL}, indicates that all estimated random effects are 
#'   conditioned on. To refrain from conditioning on any random effect parameters,
#'   specify \code{NA} or \code{~0}. The \code{newdata} argument may include new
#'   \emph{levels} of the grouping factors (for example new patient IDs),
#'   in which case the resulting posterior predictions marginalise over the 
#'   relevant random effects. See the \strong{Details} section below.
#' @param fun An optional function to apply to the predicted longitudinal response.
#'  \code{fun} is found by a call to \code{\link{match.fun}} and so can be 
#'   specified as a function object, a string naming a function, etc.
#' @param return_matrix A logical. If \code{TRUE} then a \code{draws} by 
#'   \code{nrow(newdata)} matrix is returned which contains all the actual
#'   simulations or draws from the posterior predictive distribution. Otherwise
#'   if \code{return_matrix} is set to \code{FALSE} (the default) then a 
#'   data frame is returned, as described in the \strong{Value} section below.
#' @param draws An integer indicating the number of draws to to use in
#'   generating the posterior predictions. The default
#'   and maximum number of draws is the size of the posterior sample.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param offset Not currently used. A vector of offsets. 
#'   Would only be required if \code{newdata} is specified and an \code{offset}  
#'   argument was specified when fitting the model, but offsets are not currently
#'   implemented for \code{stan_jm}.
#' @param ... Currently unused.
#' 
#' @details 
#'   If the user wishes to obtain predictions using the \code{newdata} argument
#'   then several things need to be considered: \cr 
#'   \cr
#'   First, if you wish to obtain longitudinal predictions for "new" individuals, 
#'   meaning those who were \strong{not} part of the dataset used to estimate
#'   the model, then you will likely want to marginalise over the distribution 
#'   of the individual-level random effects.
#'   To ensure that this happens, you must ensure that the IDs provided in the 
#'   \code{id_var} column of \code{newdata} do \strong{not} coincide
#'   with the IDs of individuals who were used in estimating the model. 
#'   Otherwise the predictions will be obtained using draws of the random effects 
#'   for the specific individual in the estimation data with the matching ID. 
#'   (Note that in the situation where you do want to obtain predictions for a given
#'   individual who was used in the estimation data but using new values for their 
#'   covariates, for example changing their treatment code or predicting at times 
#'   other than their actual observation times, then you could do this by specifying   
#'   the relevant individuals ID in the \code{id_var} column of \code{newdata}). \cr
#'   \cr
#'   Second, if the individuals in \code{newdata} were known to have survived beyond the 
#'   latest time provided for them in the \code{time_var} column of \code{newdata},  
#'   then this information can be provided by including a column in \code{newdata} 
#'   named "last_time". The values provided in the "last_time" column
#'   of \code{newdata} should be constant within individuals. The value provided
#'   will only be used to calculate the range for interpolation or extrapolation
#'   of time (if \code{interpolate = TRUE} or \code{extrapolate = TRUE}) as well 
#'   as being used in subsequent plots of the longitudinal trajectorie(s) to 
#'   indicate the last known survival time for the individual (assuming it was
#'   later than the time of the last predicted longitudinal response). \cr
#'   \cr
#'   Lastly, if any variables were transformed (e.g. rescaled) in the data 
#'   used to fit the model, then these variables must also be transformed in 
#'   \code{newdata}. This only applies if variables were transformed before 
#'   passing the data to one of the modeling functions and \emph{not} if 
#'   transformations were specified inside the model formula. Also see the 
#'   \strong{Note} section below for a note about using the \code{newdata} 
#'   argument with binomial models.
#' 
#' @return If \code{return_matrix = FALSE}, then a data frame of class  
#'   \code{predict.stanjm}. The data frame includes a column for the median of the
#'   posterior predictions of the mean longitudinal response (\code{yfit}),
#'   a column for each of the lower and upper limits of the uncertainty interval
#'   corresponding to the posterior predictions of the mean longitudinal response 
#'   (\code{ci_lb} and \code{ci_ub}), and a column for each of the lower and upper
#'   limits of the prediction interval corresponding to the posterior predictions
#'   of the (raw) longitudinal response. The data frame also includes columns for
#'   the subject ID variable, and each of the predictor variables. The returned
#'   object also includes a number of attributes. \cr
#'   \cr
#'   If \code{return_matrix = TRUE}, then a \code{draws} by \code{nrow(newdata)} 
#'   matrix. Each row of the matrix is a vector of posterior predictions of the 
#'   (raw) longitudinal response, generated using a single draw of the model 
#'   parameters from the posterior distribution.
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
#'   the number of trials. For more guidance, see the \pkg{rstanarm} package's
#'   \emph{How to Use the rstanarm Package} vignette.
#' 
#' @seealso \code{\link{plot.predict.stanjm}} for plotting the estimated longitudinal 
#'   trajectories, \code{\link{pp_check}} for graphical posterior predictive checks
#'   of the longitudinal submodel(s), as well as \code{\link{posterior_survfit}}
#'   for generating posterior predictions for the event submodel.
#'   
#' @examples
#' \donttest{
#'   # Run example model if not already loaded
#'   if (!exists("examplejm")) example(examplejm)
#'   
#'   # Obtain subject-specific predictions for all individuals 
#'   # in the estimation dataset
#'   preddat1 <- posterior_predict(examplejm)
#'   head(preddat1)
#'   
#'   # Obtain subject-specific predictions only for a few selected individuals
#'   preddat2 <- posterior_predict(examplejm, ids = c(1,3,8))
#'   
#'   # If we wanted to obtain subject-specific predictions in order to plot the 
#'   # longitudinal trajectories, then we might want to ensure a full trajectory 
#'   # is obtained by interpolating and extrapolating time. We can then use the 
#'   # generic plot function to plot the subject-specific predicted trajectories
#'   # for the first three individuals.
#'   preddat3 <- posterior_predict(examplejm, interpolate = TRUE, extrapolate = TRUE)
#'   head(preddat3) # predictions at additional time points compared with preddat1 
#'   plot(preddat3, ids = 1:3)
#'   
#'   # If we wanted to extrapolate further in time, but decrease the number of 
#'   # discrete time points at which we obtain predictions for each individual, 
#'   # then we could specify a named list in the 'control' argument
#'   preddat4 <- posterior_predict(examplejm, interpolate = TRUE, extrapolate = TRUE,
#'                 control = list(int_points = 10, ext_points = 10, ext_prop = 0.5))
#'   
#'   # Alternatively we may want to estimate the marginal longitudinal
#'   # trajectory for a given set of covariates. To do this, we can pass
#'   # the desired covariate values in a new data frame (however the only
#'   # covariate in our fitted model was the time variable, year). To make sure  
#'   # that we marginalise over the random effects, we need to specify an ID value
#'   # which does not correspond to any of the individuals who were used in the
#'   # model estimation. (The marginal prediction is obtained by generating 
#'   # subject-specific predictions using a series of random draws from the random 
#'   # effects distribution, and then integrating (ie, averaging) over these. 
#'   # Our marginal prediction will therefore capture the between-individual 
#'   # variation associated with the random effects.)
#'   nd <- data.frame(id = rep("new1", 11), year = (0:10 / 2))
#'   preddat5 <- posterior_predict(examplejm, newdata = nd)
#'   head(preddat5)  # note the greater width of the uncertainty interval compared 
#'                   # with the subject-specific predictions in preddat1, preddat2, etc
#'   
#'   # Alternatively, we could have estimated the "marginal" trajectory by 
#'   # ignoring the random effects (ie, assuming the random effects were set 
#'   # to zero). This will generate a predicted longitudinal trajectory only
#'   # based on the fixed effect component of the model. In essence, for a 
#'   # linear mixed effects model (ie, a model that uses an identity link 
#'   # function), we should obtain a similar point estimate ("yfit") to the
#'   # estimates obtained in preddat5 (since the mean of the estimated random effects
#'   # distribution will be approximately 0). However, it is important to note that
#'   # the uncertainty interval will be much more narrow, since it completely
#'   # ignores the between-individual variability captured by the random effects.
#'   # Further, if the model uses a non-identity link function, then the point
#'   # estimate ("yfit") obtained only using the fixed effect component of the
#'   # model will actually provide a biased estimate of the marginal prediction.
#'   # Nonetheless, to demonstrate how we can obtain the predictions only using 
#'   # the fixed effect component of the model, we simply specify 're.form = NA'. 
#'   # (We will use the same covariate values as used in the prediction for 
#'   # example for preddat5).
#'   preddat6 <- posterior_predict(examplejm, newdata = nd, re.form = NA)
#'   head(preddat6)  # note the much narrower ci, compared with preddat5
#' }
#' 
posterior_predict <- function(object, m = 1, newdata = NULL, 
                              interpolate = FALSE, extrapolate = FALSE,
                              prob = 0.95, ids, control,
                              re.form = NULL, fun = NULL, return_matrix = FALSE,
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

  # obsdata: observed data to return to user
  # newX: design matrix used for predictions
  if (is.null(newdata)) {      
    obsdata <- newX <- as.data.frame(model.frame(object)[[m]])
  } else {
    obsdata <- newX <- as.data.frame(newdata)
  }  
  
  # Error checks
  if (any(is.na(newX))) 
    stop("Currently NAs are not allowed in 'newdata'.")
  if (!id_var %in% names(newX))
    stop("The ID variable (", id_var, ") from the original fitted model ",
         "must appear in 'newdata'. To predict for new individuals ",
         "you should enter values for the ID variable that do not coincide ",
         "with those individuals used to fit the model. See the help file ",
         "for more details.", call. = FALSE)
  if (!time_var %in% names(newX))
    stop("The time variable from the original fitted model must appear ",
         "in 'newdata'.", call. = FALSE) 
  
  # Find last known survival time for each individual
  if (is.null(newdata)) { 
    # user did not provide newdata
    last_time <- object$eventtime
  } else {
    if ("last_time" %in% names(newX)) { 
      # user provided newdata with last_time column
      last_time <- tapply(newX[["last_time"]], newX[[id_var]], FUN = max)
    } else { 
      # user provided newdata but no last_time column, last_time inferred from time_var
      last_time <- tapply(newX[[time_var]], newX[[id_var]], FUN = max)
    }
  }
  
  # If user specified to predict only for a subset of IDs
  if (!is.null(ids)) {
    newX      <- newX[newX[[id_var]] %in% ids, ]
    obsdata   <- obsdata[obsdata[[id_var]] %in% ids, ]
    last_time <- last_time[as.character(ids)]
  }
  id_list <- names(last_time)
  if ((!is.null(newdata)) &&
    any(id_list %in% unique(model.frame(object)[[m]][[id_var]])))
      warning("Some of the IDs in the 'newdata' correspond to individuals in the ",
              "estimation dataset. Please be sure you want to obtain subject-",
              "specific predictions using the estimated random effects for those ",
              "individuals. If you instead meant to marginalise over the distribution ",
              "of the random effects, then please make sure the ID values do not ",
              "correspond to individuals in the estimation dataset.", immediate. = TRUE)
  
  # User specified interpolation or extrapolation
  if (interpolate || extrapolate) {  
    if (return_matrix) 
      stop("'return_matrix' cannot be TRUE if either 'interpolate' or ",
           "'extrapolate' is set to TRUE.")
    control_defaults <- list(int_points = 15, ext_points = 15, 
                             ext_distance = NULL, ext_prop = 0.2)
    if (missing(control)) {
      control <- control_defaults 
    } else if (!is.list(control)) {
      stop("'control' should be a named list.")
    } else {  # user specified control list
      if (!length(control)) {
        control <- control_defaults  
      } else {
        nms <- names(control)
        allowed_nms <- c("int_points", "ext_points", "ext_distance", "ext_prop")
        if (any(!nms %in% allowed_nms))
          stop(paste0("'control' list can only contain the following named arguments: ",
                      paste(allowed_nms, collapse = ", ")))
        if (all(c("ext_distance", "ext_prop") %in% nms))
          stop("'control' list cannot include both 'ext_distance' and 'ext_prop'.")
        if (is.null(control$int_points)) 
          control$int_points <- control_defaults$int_points   
        if (is.null(control$ext_points)) 
          control$ext_points <- control_defaults$ext_points  
        if (is.null(control$ext_distance) && is.null(control$ext_distance)) 
          control$ext_prop <-  control_defaults$ext_prop  
      }
    }
    prop <- control$ext_prop
    iinc <- control$int_points
    einc <- control$ext_points    

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
        control$ext_distance
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
    class(time_seq[[id_var]]) <- class(newX[[id_var]])
    newX[[time_var]] <- as.numeric(newX[[time_var]]) # ensures no rounding during data.table merge
    newX <- data.table::data.table(newX, key = c(id_var, time_var))
    newX <- newX[data.table::SJ(time_seq[[id_var]], time_seq[[time_var]]),
                 roll = TRUE, rollends = c(TRUE, TRUE)]  
  } else { # no interpolation or extrapolation
    if (missing(control)) control <- NULL
  }

  dat <-
    pp_data(object,
            m = m,
            newdata = newX,
            re.form = re.form,
            offset = offset,
            ...)
  ppargs <- pp_args(object, data = pp_eta(object, dat, m, draws), m)
  if (is.binomial(family(object)[[m]]$family))
    ppargs$trials <- pp_binomial_trials(object, newdata)

  ppfun <- pp_fun(object, m)
  ytilde <- do.call(ppfun, ppargs)
  if (!is.null(newX) && nrow(newX) == 1L) 
    ytilde <- t(ytilde)
  if (!is.null(fun)) 
    ytilde <- do.call(fun, list(ytilde))
  if (return_matrix) 
    return(ytilde)  # return S * N matrix, instead of data frame
  
  ytilde_med <- apply(ytilde, 2, median)
  ytilde_lb <- apply(ytilde, 2, quantile, (1 - prob)/2) 
  ytilde_ub <- apply(ytilde, 2, quantile, (1 + prob)/2) 

  mutilde <- ppargs$mu
  if (!is.null(newX) && nrow(newX) == 1L) 
    mutilde <- t(mutilde)
  mutilde_med <- apply(mutilde, 2, median)
  mutilde_lb <- apply(mutilde, 2, quantile, (1 - prob)/2) 
  mutilde_ub <- apply(mutilde, 2, quantile, (1 + prob)/2)   
    
  Terms <- terms(model.frame(object)[[m]])
  yvar <- rownames(attr(Terms, "factors"))[attr(Terms, "response")]
  xvars <- rownames(attr(Terms, "factors"))[-attr(Terms, "response")]
  dat <- if (length(xvars)) as.data.frame(newX)[, xvars, drop = FALSE] else NULL
  dat[[id_var]] <- factor(dat[[id_var]])
  out <- data.frame(cbind(dat, yfit = mutilde_med, 
                          ci_lb = mutilde_lb, ci_ub = mutilde_ub,
                          pi_lb = ytilde_lb, pi_ub = ytilde_ub))
  
  class(out) <- c("predict.stanjm", "data.frame")
  structure(out, observed_data = obsdata, last_time = last_time,
            y_var = yvar, id_var = id_var, time_var = time_var,
            interpolate = interpolate, extrapolate = extrapolate, 
            control = control, draws = draws, fun = fun, seed = seed)  
}


# call to rstanarm, to get the function to draw from the 
# appropriate posterior predictive distribution for marker m
#
# @param object A stanjm object
# @param m Integer specifying the longitudinal submodel
pp_fun <- function(object, m) {
  suffix <- family(object)[[m]]$family
  utils::getFromNamespace(paste0(".pp_", suffix), "rstanarm")
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
  if (is.gaussian(famname)) {
    args$sigma <- stanmat[, paste0("Long", m, "|sigma")]
  } else if (is.gamma(famname)) {
    args$shape <- stanmat[, paste0("Long", m, "|shape")]
  } else if (is.ig(famname)) {
    args$lambda <- stanmat[, paste0("Long", m, "|lambda")]
  } else if (is.nb(famname)) {
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
  S <- posterior_sample_size(object)
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
    eta <- linear_predictor.matrix(beta, x, data$offset)
  } else {
    stanmat <- as.matrix(object$stanfit)
    nms <- collect_nms(colnames(stanmat), M)   
    beta <- stanmat[, nms$y[[m]], drop = FALSE]
    if (some_draws) 
      beta <- beta[samp, , drop = FALSE]
    eta <- linear_predictor.matrix(beta, x, data$offset)
    b <- stanmat[, nms$y_b[[m]], drop = FALSE]
    if (some_draws) 
      b <- b[samp, , drop = FALSE]
    if (is.null(data$Z_names)) {
      b <- b[, !grepl("_NEW_", colnames(b), fixed = TRUE), drop = FALSE]
    } else {
      b <- pp_b_ord(b, data$Z_names)
    }
    eta <- eta + as.matrix(b %*% data$Zt)
  }
  nlist(eta, stanmat)
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
