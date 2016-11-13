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

#' Estimate marginal or subject-specific survival probabilities
#' 
#' The posterior predictive distribution is the distribution of the outcome 
#' implied by the model after using the observed data to update our beliefs 
#' about the unknown parameters in the model. This function allows us to 
#' generate estimated survival probabilities (either subject-specific, or
#' by marginalising over the distribution of the random effects) based on
#' draws from the posterior predictive distribution. In both the "subject-specifc"
#' and "marginal" situations, the predicted survival probabilities will still be 
#' \emph{conditional} on observed values of the fixed effect covariates 
#' in the longitudinal and event submodels (ie, the predictions will be obtained 
#' using either the design matrices used in the original \code{\link{stan_jm}} model
#' call, or using the covariate values provided in the \code{newdata} argument). However, 
#' if you wish to also average over the observed distribution of the fixed effect 
#' covariates then this is also possible -- however we refer to these
#' as standardised survival probabilties -- see the \code{standardise} 
#' argument below.
#' 
#' @export
#' @templateVar stanjmArg object
#' @template args-stanjm-object
#' 
#' @param newdata Optionally, a new data frame in which to look 
#'   for variables with which to predict. If omitted, the model matrices are used. 
#'   If new data is provided, then it should contain the covariate values needed
#'   for all longitudinal submodel(s) and the event submodel. There is only
#'   allowed to be one row of data for each individual in \code{newdata}, that
#'   is, time-varying covariates are not allowed in the prediction dataset. Also 
#'   note that if \code{newdata} is provided, then the \code{times} argument
#'   must also be specified. See the \strong{Details} section for further 
#'   important details that need to be considered when specifying \code{newdata}.
#' @param extrapolate A logical specifying whether to extrapolate the estimated 
#'   survival probabilities beyond the times specified in the \code{times} argument.
#'   If \code{TRUE} then the extrapolation can be further controlled using
#'   the \code{control} argument.
#' @param control A named list with parameters controlling extrapolation 
#'   of the estimated survival function when \code{extrapolate = TRUE}. The list
#'   can contain one or more of the following named elements: \cr
#'   \describe{
#'     \item{\code{ext_points}}{a positive integer specifying the number of  
#'     discrete time points at which to calculate the forecasted survival 
#'     probabilities. The default is 15.}
#'     \item{\code{ext_prop}}{a positive scalar between 0 and 1 specifying the 
#'     amount of time across which to forecast the estimated survival function,
#'     represented as a proportion of the total observed follow up time for each
#'     individual. For example specifying \code{ext_prop = 0.2} means that for an
#'     individual for whom the latest of their measurement, event or censoring times
#'     was 10 years, their estimated survival function will be extrapolated 
#'     out to 12 years (i.e. 10 + (0.2 * 10)). The default value is 0.2.}
#'     \item{\code{ext_distance}}{a positive scalar specifying the amount of time 
#'     across which to forecast the estimated survival function, represented 
#'     in units of the time variable \code{time_var} (from fitting the model). 
#'     This cannot be specified if \code{ext_prop} is specified.}
#'     \item{\code{condition}}{a logical specifying whether the estimated 
#'     subject-specific survival probabilities at time \code{t} should be 
#'     conditioned on survival up to a fixed time point \code{u}. The default 
#'     is to condition on the latest observation time for each individual 
#'     (taken to be the event or censoring time if \code{newdata} is not 
#'     specified, or the value of the \code{times} if \code{newdata} is 
#'     specified but no \code{last_time} is provided in the \code{control} 
#'     list, or otherwise the times provided in the \code{last_time} element
#'     of the \code{control} list).}
#'     \item{\code{last_time}}{a scalar, a numeric vector or a character string 
#'     specifying the last known survival time for each individual for whom
#'     conditional predictions are being obtained. Should only be specified if
#'     \code{newdata} is provided. If \code{last_time} is not provided then
#'     the default is to use the value provided in the \code{times} argument.
#'     A scalar will use the same last time for each individual in \code{newdata}.
#'     A numeric vector should provide a last time for each individual in 
#'     \code{newdata}. A character string should name a column in \code{newdata}
#'     in which to look for the last times.} 
#'   }
#' @param ids An optional vector specifying a subset of IDs for whom the 
#'   predictions should be obtained. The default is to predict for all individuals
#'   who were used in estimating the model or, if \code{newdata} is specified,
#'   then all individuals contained in \code{newdata}.
#' @param prob A scalar between 0 and 1 specifying the width to use for the 
#'   uncertainty interval (sometimes called credible interval) for the predictions. 
#'   For example \code{prob = 0.95} (the default) means that the 2.5th and 97.5th  
#'   percentiles will be provided.
#' @param times A scalar, a numeric vector or a character string specifying the values of 
#'   \code{time_var} in the original model at which to obtain the estimated 
#'   survival probabilities. If \code{newdata} is not provided, then the 
#'   \code{times} argument is optional, and if not provided then \code{times} 
#'   will default to the last known event or censoring time for each individual.
#'   If \code{newdata} is provided, then the \code{times} argument must be 
#'   provided and the user has several options. If a scalar is supplied, then 
#'   the survival probabilities will all be calculated at the time specified. 
#'   If a numeric vector is supplied then it should be the same length as the 
#'   number of new individuals for whom predictions are being obtained, with 
#'   the \eqn{i}th element of the vector specifying the time at which to 
#'   calculate the survival probability for the \eqn{i}th individual in the 
#'   \code{newdata}. If a character string is supplied then it should be the 
#'   name of the column in \code{newdata} which provides the times at which to
#'   calculate the survival probabilities.
#' @param standardise A logical specifying whether the estimated 
#'   subject-specific survival probabilities should be averaged
#'   across all individuals for whom the subject-specific predictions are 
#'   being obtained. This can be used to average over the covariate profile of
#'   either the individuals used in estimating the model, or the covariate 
#'   profiles of the individuals provided in \code{newdata}. This approach of
#'   averaging across the observed distribution of the covariates is sometimes
#'   referred to as a "standardised" survival curve. If \code{TRUE}, then the 
#'   \code{times} argument must be specified and must be of length 1 or specify
#'   a column in \code{newdata} that contains times which are constant across 
#'   individuals.
#' @param draws An integer indicating the number of MCMC draws to return. The default
#'   and maximum number of draws is the size of the posterior sample.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param offset Not currently used. A vector of offsets. Would  
#'   only be required if \code{newdata} is specified and an \code{offset}  
#'   argument was specified when fitting the model, but offsets are not currently
#'   implemented for \code{stan_jm}.
#'
#' @details 
#'   If the user wishes to obtain predictions using the \code{newdata}
#'   argument then several things need to be considered: \cr 
#'   \cr
#'   First, if you wish to obtain survival probabilities for "new" individuals, 
#'   meaning those who were \strong{not} part of the dataset used to estimate
#'   the model, then you will likely want to marginalise over the distribution 
#'   of the individual-level random effects.
#'   To ensure that this happens, you must ensure that the IDs provided in the 
#'   \code{id_var} column of \code{newdata} do 
#'   \strong{not} coincide with the IDs of individuals who were used in estimating  
#'   the model. Otherwise the predictions will be obtained using draws of the random  
#'   effects for the specific individual in the estimation data with the matching ID. 
#'   (Note that in the situation where you do want to obtain predictions for a given
#'   individual who was used in the estimation data but using new values for their 
#'   covariates, for example changing their treatment code or predicting at times 
#'   other than their actual observation times, then you could do this by specifying   
#'   the relevant individuals ID in the \code{id_var} columns of \code{newdata}. \cr
#'   \cr
#'   Second, if any variables were transformed (e.g. rescaled) in the data 
#'   used to fit the model, then these variables must also be transformed in 
#'   \code{newdata}. This only applies if variables  
#'   were transformed before passing the data to one of the modeling functions and  
#'   \emph{not} if transformations were specified inside the model formula. Also  
#'   see the \strong{Note} section in \code{\link{posterior_predict}} for a note  
#'   about using the \code{newdata} argument with binomial models.
#'    
#' @return A data frame of class \code{survfit.stanjm}. The data frame includes 
#'   columns for each of the following: \cr
#'   \describe{
#'     \item{the median of the posterior predictions of the estimated survival
#'     probabilities (\code{survfit})}
#'     \item{each of the lower and upper limits of the corresponding uncertainty 
#'     interval for the estimated survival probabilities (\code{ci_lb} and 
#'     \code{ci_ub})}
#'     \item{a subject identifier (\code{id_var}), unless standardised survival
#'     probabilities were estimated
#'     \item{the time that the estimated survival probability is calculated for 
#'     (\code{time_var})
#'     \item{if \code{newdata} was provided, then the values of each of the
#'     covariates used in the calculation of the estimated survival probabilities}
#'   }
#'   The returned object also includes a number of attributes.
#' 
#' @seealso \code{\link{ps_check}} for for graphical checks of the estimated 
#'   survival function, and \code{\link{posterior_predict}} for estimating the
#'   marginal or subject-specific longitudinal trajectories.
#'   
#' @examples
#' \donttest{
#'   # Run example model if not already loaded
#'   if (!exists("examplejm")) example(examplejm)
#'   
#'   # Obtain subject-specific survival probabilities for a few
#'   # selected individuals in the estimation dataset who were  
#'   # known to survive up until their censoring time 
#'   head(pbcSurv_subset[pbcSurv_subset$status == 0,])
#'   ps1 <- posterior_survfit(examplejm, ids = c(7,13,16))
#'   head(ps1)
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
#' }
#'  
posterior_survfit <- function(object, newdata = NULL,  
                              extrapolate = TRUE, control = list(), 
                              prob = 0.95, ids,
                              times = NULL, standardise = FALSE, 
                              draws = NULL, seed = NULL, offset = NULL, ...) {
  validate_stanjm_object(object)
  M <- object$n_markers
  id_var <- object$id_var
  time_var <- object$time_var
  if (missing(ids)) ids <- NULL
  if (standardise) {
    if (is.null(times) || (length(times) > 1L))
      stop("'times' must be specified and be of length 1 in order to obtain ",
           "standardised survival probabilities. (The subject-specific survival ",
           "probabilities will be calculated at the specified time point, and ",
           "then averaged).")
    if (condition)
      stop("'condition' cannot be set to TRUE if standardised survival ",
           "probabilities are requested.")
  }
  if (!is.null(seed)) 
    set.seed(seed)
  
  # Construct ndE and ndL
  if (is.null(newdata)) {
    ndL <- model.frame(object)[1:M]
    ndE <- model.frame(object)$Event
    id_list <- unique(ndE[[id_var]])
    if (is.null(times)) 
      times <- object$eventtime
    if (!is.null(control$last_time)) {
      stop("'last_time' should not be provided when 'newdata' is NULL, since ",
           "times are taken to be the event or censoring time for each individual")
    } else {
      last_time <- object$eventtime
    }
  } else {  # user specified newdata
    if (any(is.na(newdata))) 
      stop("Currently NAs are not allowed in 'newdata'.")
    if (!id_var %in% colnames(newdata))
      stop("id_var from the original model call must appear 'newdata'.")
    ndL <- lapply(1:M, function(m) as.data.frame(newdata))
    ndE <- as.data.frame(newdata)
    id_list <- unique(ndE[[id_var]])
    if (!identical(length(id_list), nrow(ndE)))
      stop("'newdata' should only contain one row per individual, since ",
           "time varying covariates are not allowed in the prediction data.")
    if (any(id_list %in% unique(model.frame(object)[[1]][[id_var]])))
      warning("Some of the IDs in 'newdata' correspond to individuals in the ",
              "estimation dataset. Please be sure you want to obtain subject-",
              "specific predictions using the estimated random effects for those ",
              "individuals. If you instead meant to marginalise over the distribution ",
              "of the random effects, then please make sure the ID values do not ",
              "correspond to individuals in the estimation dataset.", immediate. = TRUE)
    if (is.null(times)) 
      stop("'times' cannot be NULL if newdata is specified.")
    if (!is.vector(times) || !is.numeric(times))
      stop("'times' should be a numeric vector.")
    if (is.null(control$last_time)) last_time <- NULL
  }

  # If user specified to predict only for a subset of IDs
  if (!is.null(ids)) {
    if (!all(ids %in% id_list))
      stop("Some 'ids' do not appear in the data.")
    ndE <- ndE[ndE[[id_var]] %in% ids,]
    ndL <- lapply(ndL, function(x) x[x[[id_var]] %in% ids, ])
    id_list <- id_list[id_list %in% ids]
    if (is.null(newdata)) {
      times <- times[as.character(ids)]
      last_time <- last_time[as.character(ids)]
    }
  }
  if (length(times) == 1L)
    times <- rep(times, length(id_list))
  if (is.null(last_time)) 
    last_time <- times  
  if (!identical(length(times), length(id_list)))
    stop(paste0("'times' vector should be of length 1 or length equal to the ",
                "number of individuals for whom predictions are being obtained (",
                length(id_list), ")."))  
  
  # User specified extrapolation
  if (extrapolate) {
    control_defaults <- list(ext_points = 15, ext_distance = NULL, 
                             ext_prop = 0.2, condition = TRUE, last_time = NULL) 
    if (!length(control)) {
      control <- control_defaults 
    } else if (!is.list(control)) {
      stop("'control' should be a named list.")
    } else {  # user specified control list
      if (!length(control)) {
        control <- control_defaults  
      } else {
        nms <- names(control)
        allowed_nms <- c("ext_points", "ext_distance", "ext_prop", "condition", "last_time")
        if (any(!nms %in% allowed_nms))
          stop(paste0("'control' list can only contain the following named arguments: ",
                      paste(allowed_nms, collapse = ", ")))
        if (all(c("ext_distance", "ext_prop") %in% nms))
          stop("'control' list cannot include both 'ext_distance' and 'ext_prop'.")
        if (is.null(control$ext_points)) 
          control$ext_points <- control_defaults$ext_points  
        if (is.null(control$ext_distance) && is.null(control$ext_distance)) 
          control$ext_prop <-  control_defaults$ext_prop
        if (is.null(control$condition))
          control$condition <- TRUE
      }
    }
    prop <- control$ext_prop
    inc <- control$ext_points 
    # Note: prop assumes all individuals entered at time 0
    dist <- if (!is.null(prop)) prop * (last_time - 0) else
      control$ext_distance
    time_seq <- lapply(0:inc, function(x, t) t + dist * (x / inc), 
                       t = times)    
  } else { # no extrapolation
    if (missing(control)) control <- NULL
    time_seq <- list(times)
  }    

  surv <- lapply(time_seq, function(x) {
    dat <- ps_data(object,
                  newdataEvent = ndE,
                  newdataLong = ndL,
                  ids = id_list,
                  t = x, 
                  id_var = id_var,
                  time_var = time_var,
                  ...) 
    surv <- ps_survcalc(object, dat, draws)
    if (!is.null(newdata) && nrow(newdata) == 1L) 
      surv <- t(surv)
    # set survprob matrix at time 0 to S(t) = 1 
    # (otherwise some NaN possible due to numerical inaccuracies)
    surv[,(x == 0)] <- 1
    if (standardise) {
      surv <- matrix(rowMeans(surv), ncol = 1)
      dimnames(surv) <- list(iterations = NULL, "standardised_survprob")
    } else {
      dimnames(surv) <- list(iterations = NULL, ids = id_list)
    }
    surv
  })
  # Optionally condition on survprob matrix at last_time
  if (extrapolate && control$condition) {
    cond_dat <- ps_data(object,
                   newdataEvent = ndE,
                   newdataLong = ndL,
                   ids = id_list,
                   t = last_time, 
                   id_var = id_var,
                   time_var = time_var,
                   ...) 
    cond_surv <- ps_survcalc(object, cond_dat, draws)
    if (!is.null(newdata) && nrow(newdata) == 1L) 
      cond_surv <- t(cond_surv)
    # set survprob matrix at time 0 to S(t) = 1 
    # (otherwise some NaN possible due to numerical inaccuracies)
    cond_surv[,(last_time == 0)] <- 1    
    surv <- lapply(surv, function(x) {
      vec <- x / cond_surv
      vec[is.na(vec)] <- 1
      vec})
  }
  
  # Summarise posterior draws to get median and ci
  out <- do.call("rbind", 
    lapply(seq_along(surv), function(x, limits, id_list, standardise, 
                                     id_var, time_var) {
      surv_med <- apply(surv[[x]], 2, median)
      surv_lb <- apply(surv[[x]], 2, quantile, (1 - prob)/2) 
      surv_ub <- apply(surv[[x]], 2, quantile, (1 + prob)/2)  
      out <- cbind(IDVAR = if (!standardise) id_list, 
                   TIMEVAR = if (!standardise) time_seq[[x]] else unique(time_seq[[x]]),
                   surv_med, surv_lb, surv_ub)
      out
    }, limits = limits, id_list = id_list, standardise = standardise, 
       id_var = id_var, time_var = time_var))
  rownames(out) <- NULL
  colnames(out) <- c(if ("IDVAR" %in% colnames(out)) id_var,
                     time_var, "survpred", "ci_lb", "ci_ub")    
  if (id_var %in% colnames(out)) {  # data has id column -- sort by id and time
    out <- out[order(out[, id_var, drop = F], out[, time_var, drop = F]), , drop = F]
  } else { # data does not have id column -- sort by time only
    out <- out[order(out[, time_var, drop = F]), , drop = F]
  }
  out <- data.frame(out)
  class(out) <- c("survfit.stanjm", "data.frame")
  structure(out, id_var = id_var, time_var = time_var,
            extrapolate = extrapolate, control = control, 
            standardise = standardise, ids = id_list, 
            draws = draws, seed = seed, offset = offset)
}


# create matrix of posterior survival probabilities at time t
# 
# @param object stanjm object
# @param data output from ps_data()
# @param draws number of draws
# @return matrix of survival probabilities at time t, with the S
#   rows corresponding to different MCMC draws of the parameters 
#   from the posterior, and each column corresponding to Npat
#   individuals in the new data
ps_survcalc <- function(object, data, draws = NULL) {
  M <- object$n_markers
  e_xQ <- data$e_xQ
  y_xQ <- data$y_xQ
  y_zQ <- data$y_zQ
  t    <- data$t
  tQ   <- data$tQ
  t_and_tQ  <- c(list(t), tQ)
  Npat      <- length(t)
  quadnodes <- length(tQ)
  Q         <- quadnodes + 1
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
  stanmat <- as.matrix(object$stanfit)
  nms <- collect_nms(colnames(stanmat), M)
  
  # Longitudinal submodels
  y_beta <- lapply(1:M, function(m) {
    mat <- stanmat[, nms$y[[m]], drop = FALSE]
    if (some_draws) 
      mat <- mat[samp, , drop = FALSE]
    mat
  }) 
  eta_long <- lapply(seq(M), function(m) {
    lapply(seq(Q), function(q)
      linear_predictor.matrix(
        y_beta[[m]], y_xQ[[m]][[q]], data$offset)) 
  }) 
    
  y_b <- lapply(seq(M), function(m) {
    mat <- stanmat[, nms$y_b[[m]], drop = FALSE]
    if (some_draws) 
      mat <- mat[samp, , drop = FALSE]
    pp_b_ord(mat, y_zQ[[m]]$Z_names)
  })
  eta_long <- lapply(seq(M), function(m) {
    lapply(seq(Q), function(q)
      eta_long[[m]][[q]] + as.matrix(y_b[[m]] %*% y_zQ[[m]]$Zt[[q]]))
  })

  # Event submodel
  e_beta <- stanmat[, nms$e, drop = FALSE]
  if (some_draws) 
    e_beta <- e_beta[samp, , drop = FALSE] 
  eta_event <- lapply(seq(Q), function(q)
    linear_predictor.matrix(e_beta, e_xQ[[q]], offset = NULL))
  
  # Association structure
  assoc <- object$assoc
  if (any(unlist(assoc))) {
    a_beta <- stanmat[, nms$a, drop = FALSE]
    if (some_draws) 
      a_beta <- a_beta[samp, , drop = FALSE]
    mark <- 1
    for (m in 1:M) {
      if (assoc$etavalue[m]) {
        eta_event <- lapply(seq(Q), function(q)
          eta_event[[q]] + a_beta[, mark] * eta_long[[m]][[q]]) 
        mark <- mark + 1  
      } 
      if (assoc$etaslope[m]) {
        #etaslope_long[[m]] <- NULL  # !!! needs calculation of derivative
        #eta_event <- eta_event + a_beta[, mark] * etaslope_long[[m]] 
        mark <- mark + 1  
      }
      if (assoc$muvalue[m]) {
        invlink <- family(object)[[m]]$linkinv
        eta_event <- lapply(seq(Q), function(q)
          eta_event[[q]] + a_beta[, mark] * invlink(eta_long[[m]][[q]])) 
        mark <- mark + 1  
      }
      if (assoc$muslope[m]) {
        #muslope_long[[m]] <- NULL  # !!! needs calculation of derivative
        #eta_event <- eta_event + a_beta[, mark] * muslope_long[[m]] 
        mark <- mark + 1  
      }    
    }
    if (any(assoc$shared_b)) {
      # !!!
    }    
  }

  # Baseline hazard
  if (object$base_haz$type == "weibull") {
    shape <- stanmat[, nms$e_extra, drop = FALSE]
    if (some_draws) 
      shape <- shape[samp, , drop = FALSE] 
    log_basehaz <- lapply(t_and_tQ, function(x) {
      # returns S x Npat matrix
      as.vector(log(shape)) + (shape - 1) %*% matrix(log(x), nrow = 1)
    })
  } else if (object$base_haz$type == "splines") {
    coefs <- stanmat[, nms$e_extra, drop = FALSE]
    if (some_draws) 
      coefs <- coefs[samp, , drop = FALSE] 
    log_basehaz <- lapply(t_and_tQ, function(x) {
      # returns S x Npat matrix
      coefs %*% t(splines::ns(x, knots = object$base_haz$splines_attr$knots,
                     intercept = object$base_haz$splines_attr$intercept,
                     Boundary.knots = object$base_haz$splines_attr$Boundary.knots))
    })
  }
  log_haz <- mapply(function(x,y) x + y, 
                    log_basehaz, eta_event, SIMPLIFY = FALSE)
  
  log_haz_t <- log_haz[[1]]
  log_haz_Q <- log_haz[2:Q]
  qw <- get_quadpoints(quadnodes)$weights
  qw_times_half_t <- lapply(seq(quadnodes), function(x) qw[x] * (t / 2))
  haz_Q <- lapply(log_haz_Q, exp)
  weighted_haz_Q <- lapply(seq(quadnodes), function(q) {
    # returns S x Npat matrix
    t(apply(haz_Q[[q]], 1, function(row) row * qw_times_half_t[[q]]))    
  })
  sum_weighted_haz_Q <- Reduce('+', weighted_haz_Q)
  surv_t <- exp(-sum_weighted_haz_Q)
  return(surv_t) # returns S x Npat matrix of survival probabilities at t
}








