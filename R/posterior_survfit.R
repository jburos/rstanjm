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

#' Draw from posterior predictive distribution for the event submodel
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
#' @template args-stanjm-object
#' @param newdata Optionally, a data frame in which to look for variables with 
#'   which to predict. If omitted, the model matrix is used. If \code{newdata} 
#'   is provided and any variables were transformed (e.g. rescaled) in the data 
#'   used to fit the model, then these variables must also be transformed in 
#'   \code{newdata}. This only applies if variables were transformed before 
#'   passing the data to one of the modeling functions and \emph{not} if 
#'   transformations were specified inside the model formula. Also see the Note
#'   section below for a note about using the \code{newdata} argument with with
#'   binomial models.
#' @param ids An optional character vector containing the IDs of individuals for 
#'   whom the estimated survival probabilites should be obtained. This defaults 
#'   to \code{NULL} which returns predictions for all individuals in the original 
#'   model if \code{newdataEvent} is \code{NULL}, or otherwise all individuals 
#'   in \code{newdataEvent}.
#' @param times An optional scalar or numeric vector specifying the values of 
#'   \code{time_var} in the original model at which to obtain the estimated 
#'   survival probabilities. If not specified, and \code{marginalised = FALSE}, 
#'   then the default behaviour is to use the latest observation time for each 
#'   individual (taken to be the event or censoring time if \code{newdata} is 
#'   not specified, or the lastest value of \code{time_var} if \code{newdata} 
#'   is specified) and extrapolate forward from there. If not specified, and
#'   \code{marginalised = TRUE}, then the default behaviour is to calculate
#'   the marginal survival probability at \code{n_increments} time points between
#'   baseline and the latest observation or event/censoring time for all individuals.
#' @param condition A logical specifying whether the estimated subject-specific
#'   survival probabilities at time \code{t} should be conditioned on 
#'   survival up to a fixed time point \code{u}. The default is to condition
#'   on the latest observation time for each individual (taken to be the event 
#'   or censoring time if \code{newdata} is not specified, or the lastest value of 
#'   \code{time_var} if \code{newdata} is specified). It cannot
#'   be set to \code{TRUE} if \code{marginalised} is set to \code{TRUE}.
#' @param extrapolate A logical specifying whether to extrapolate the estimated 
#'   survival function beyond the times specified in the \code{times} argument
#'   (which will default to the latest observation time -- measurement, event or
#'   censoring time -- if left equal to \code{NULL} and \code{marginalised = FALSE}).
#'   If set to \code{TRUE} then the predictions can be further controlled using
#'   the \code{extrapolate_args} argument.
#' @param extrapolate_args A named list with parameters controlling extrapolation 
#'   of the estimated survival function when \code{extrapolate = TRUE}. The list
#'   can contain the following named elements: 
#'   \code{dist}, a positive scalar 
#'   specifying the amount of time across which to forecast the estimated survival
#'   function; 
#'   \code{prop}, a positive scalar between 0 and 1 specifying the 
#'   amount of time across which to forecast the estimated survival function but
#'   represented as a proportion of the total observed follow up time for each
#'   individual;
#'   \code{increments}, a positive integer specifying the number of increments in 
#'   time at which to calculate the forecasted survival probabilities. 
#'   See the \strong{Examples} below.
#' @param marginalised A logical specifying whether the estimated 
#'   subject-specific survival probabilities should be averaged (marginalised)
#'   across all individuals for whom the predictions were obtained. This can be
#'   used to obtain an estimate of the marginal survival probability. It cannot
#'   be set to \code{TRUE} if \code{condition} is set to \code{TRUE}.
#' @param draws An integer indicating the number of MCMC draws to return. The default
#'   and maximum number of draws is the size of the posterior sample.
#' @param fun An optional function to apply to the results. \code{fun} is found 
#'   by a call to \code{\link{match.fun}} and so can be specified as a function
#'   object, a string naming a function, etc.
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param ... Currently unused.
#' 
#' @return A named list with two elements, 'times' and 'survprobs'. The former
#'   is a list of times at which the survival probabilities are calculated. 
#'   The latter is a list containing the correponding survival probabilities.
#'   Each element of 'survprobs' contains a \code{draws} by \code{levels(ids)} 
#'   matrix of estimated survival probabilities. Each row of the matrix is a
#'   vector of predictions generated using a single draw of the model parameters
#'   from the posterior distribution.
#' 
#' @seealso \code{\link{ps_check}} for for graphical checks of the estimated 
#'   marginal survival function, and \code{\link{posterior_predict}} for 
#'   drawing from the posterior predictive distribution for the longitudinal 
#'   submodel.
#'   
#' @examples
#' The numeric value should represent the proportion of observed follow up time
#'   (difference between baseline and \code{time}), for example, 0.2 = 20% of the 
#'   observed follow up time. If set to \code{NULL} then survival probabilities
#'   are only calculated at \code{time} and no extrapolation is carried out.
#'   It is ignored if \code{marginalised} is set to \code{TRUE}.
#' 
posterior_survfit <- function(object, newdataEvent = NULL, newdataLong = NULL, ids,  
                              times = NULL, condition = TRUE, extrapolate = TRUE, 
                              extrapolate_args = list(dist = NULL, prop = 0.5, 
                                                      increments = 25), 
                              marginalised = FALSE, 
                              draws = NULL, fun = NULL, seed = NULL, ...) {
  validate_stanjm_object(object)
  M <- object$n_markers
  id_var <- object$id_var
  time_var <- object$time_var
  if (missing(ids)) ids <- NULL
  if (marginalised) {
    if (condition || extrapolate)
      cat("'marginalised' survival predictions were requested, therefore, ",
          "'condition' and 'extrapolate' will both automatically be set ",
          "to FALSE.")
      condition <- extrapolate <- FALSE
    if ((!is.null(times)) && (length(times) > 1L))
      stop("'times' can only be length 1 if 'marginalised = TRUE'.")
  }
  if (extrapolate) {
    if (!is.list(extrapolate_args))
      stop("'extrapolate_args' must be a named list.")
    if ((!is.null(extrapolate_args$dist)) && (extrapolate_args$dist < 0))
      stop("'extrapolate_args$dist' must be positive or NULL.")
    if ((!is.null(extrapolate_args$prop)) && (extrapolate_args$prop < 0))
      stop("'extrapolate_args$prop' must be positive or NULL.")
    if (is.null(extrapolate_args$increments) || (extrapolate_args$increments < 0))
      stop("'extrapolate_args$increments' must be a non-negative integer.")
    extrapolate_args$increments <- as.integer(extrapolate_args$increments)
    if ((is.null(extrapolate_args$prop)) && (is.null(extrapolate_args$dist)))
      stop("One of 'extrapolate_args$dist' or 'extrapolate_args$prop' ", 
           "must be specified as non-NULL if 'extrapolate = TRUE'.")    
  }
  if (!is.null(seed)) 
    set.seed(seed)
  if (!is.null(fun)) 
    fun <- match.fun(fun)
  if ((!is.null(newdataEvent)) && is.null(newdataLong))
      stop("newdataLong must also be specified if specifying newdataEvent.")
  if (!is.null(newdataLong)) {
    if (is.null(newdataEvent))
      stop("newdataEvent must also be specified if specifying newdataLong.")
    newdataEvent <- as.data.frame(newdataEvent)
    if (is.data.frame(newdataLong)) {
      newdataLong <- list(newdataLong)
    } else if (is(newdataLong, "list")) {
      newdataLong <- lapply(newdataLong, as.data.frame)
    } else {
      stop("newdataLong should be a data frame or a list of data frames.")
    }
    if (!identical(length(newdataLong), M)) {
      if (M == 1) 
        stop("newdataLong should be a data frame.") else
          stop(paste0("newdataLong appears to be a list of the incorrect length. ",
                      "It should be the same length as the number of longitudinal ",
                      "markers (M=", M))
    }
    lapply(c(newdataLong, list(newdataEvent)), function(x) {
      if (!id_var %in% colnames(x))
        stop("id_var from the original model call must appear in all ",
             "newdata data frames.")
      if (!time_var %in% colnames(x))
        stop("time_var from the original model call must appear in all ",
             "newdata data frames.")
    })    
    if (!is.null(ids)) {
      newdataEvent <- newdataEvent[newdataEvent[[id_var]] %in% ids,]
      newdataLong <- lapply(newdataLong, function(x) x[x[[id_var]] %in% ids, ])
    }    
    if (any(is.na(newdataEvent))) 
      stop("Currently NAs are not allowed in 'newdataEvent'.")
    lapply(newdataLong, function(x) 
      if (any(is.na(x)))
        stop("Currently NAs are not allowed in 'newdataLong'."))
    ids_and_times <- do.call(rbind, lapply(c(newdataLong, list(newdataEvent)), 
                       function(x) x[, c(id_var, time_var)]))
    # Latest known observation time for each individual
    lasttime <- tapply(ids_and_times[[time_var]], ids_and_times[[id_var]], FUN = max)
  } else {
    newdataLong <- model.frame(object)[1:M]
    newdataEvent <- model.frame(object)$Event
    # Latest known observation time for each individual
    lasttime <- object$eventtime
    if (!is.null(ids)) {
      newdataEvent <- newdataEvent[newdataEvent[[id_var]] %in% ids,]
      newdataLong <- lapply(newdataLong, function(x) x[x[[id_var]] %in% ids, ])
      lasttime <- lasttime[ids]
    }    
  }
  
  # Time sequence across which to generate the survival probabilities
  max_time <- max(object$eventtime)
  prop <- extrapolate_args$prop
  inc <- extrapolate_args$increments
  if (marginalised) {  
    # marginalised predictions of survival probability
    if (is.null(times)) {
      max_time <- rep(max_time, length(lasttime))
      time_seq <- lapply(0:inc, function(x) 
        (x * max_time) / inc)
    } else time_seq <- list(rep(times, length(lasttime)))
  } else {
    # subject-specific predictions
    if ((!is.null(times)) && (length(times) == 1L))
      times <- rep(times, length(lasttime))
    if (extrapolate) {  
      # extrapolate
      # Note: prop assumes all individuals entered at time 0
      # Note: still use whole follow up period for determining how 
      #       far to extrapolate even if times are provided?
      dist <- if (!is.null(prop)) prop * (lasttime - 0) else
        extrapolate_args$prop
      time_seq <- lapply(0:inc, function(x, t) t + dist * (x / inc), 
        t = if (!is.null(times)) times else lasttime)
    } else {
      # don't extrapolate
      time_seq <- if (!is.null(times)) list(times) else list(lasttime)
    }
  }

  # List of ordered ids
  id_list <- unique(newdataEvent[[id_var]])
  
  surv <- lapply(time_seq, function(x) {
    dat <- ps_data(object,
                  newdataEvent = newdataEvent,
                  newdataLong = newdataLong,
                  ids = id_list,
                  t = x, 
                  id_var = id_var,
                  time_var = time_var,
                  ...) 
    surv <- ps_survcalc(object, dat, draws)
    if (!is.null(newdataEvent) && nrow(newdataEvent) == 1L) 
      surv <- t(surv)
    # set survprob matrix at time 0 to S(t) = 1 
    # (otherwise some NaN possible due to numerical inaccuracies)
    surv[,(x == 0)] <- 1
    if (marginalised) {
      surv <- matrix(rowMeans(surv), ncol = 1)
      dimnames(surv) <- list(iterations = NULL, "marginal_survprob")
    } else {
      dimnames(surv) <- list(iterations = NULL, ids = id_list)
    }
    if (!is.null(fun))
      stop("'fun' not currently implemented.")
    #  surv <- do.call(fun, list(surv))
    surv
  })
  # If marginal survprob then time_seq is same for all ids
  if (marginalised) {
    time_seq <- lapply(time_seq, unique)
  }
  # Optionally condition on first survprob matrix (increment 0)
  if (condition) {
    surv <- lapply(surv, function(x) x / surv[[1]])
  }
  names(surv) <- names(time_seq) <- paste0("increment", 0:inc)
  #if (length(time_seq) == 1L) {
  #  time_seq <- unlist(time_seq)
  #  survprobs <- unlist(survprobs)    
  #}

  out <- list(times = time_seq, survprobs = surv)
  class(out) <- "survfit.stanjm"
  structure(out,
            marginalised = marginalised, condition = condition,
            extrapolate = extrapolate, extrapolate_args = extrapolate_args, 
            ids = id_list, draws = draws, fun = fun, seed = seed)
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
      rstanarm:::linear_predictor.matrix(
        y_beta[[m]], y_xQ[[m]][[q]], data$offset)) 
  }) 
    
  y_b <- lapply(seq(M), function(m) {
    mat <- stanmat[, nms$y_b[[m]], drop = FALSE]
    if (some_draws) 
      mat <- mat[samp, , drop = FALSE]
    rstanarm:::pp_b_ord(mat, y_zQ[[m]]$Z_names)
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
    rstanarm:::linear_predictor.matrix(e_beta, e_xQ[[q]], offset = NULL))
  
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
  if (object$base_haz == "weibull") {
    shape <- stanmat[, nms$e_extra, drop = FALSE]
    if (some_draws) 
      shape <- shape[samp, , drop = FALSE] 
    log_basehaz <- lapply(t_and_tQ, function(x) {
      # returns S x Npat matrix
      as.vector(log(shape)) + (shape - 1) %*% matrix(log(x), nrow = 1)
    })
  } else if (object$base_haz == "splines") {
    coefs <- stanmat[, nms$e_extra, drop = FALSE]
    if (some_draws) 
      coefs <- coefs[samp, , drop = FALSE] 
    log_basehaz <- lapply(t_and_tQ, function(x) {
      # returns S x Npat matrix
      coefs %*% t(ns(x, object$df)) # !!! needs to accept df or knots
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








