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

#' Bayesian joint longitudinal and time-to-event models via Stan
#' 
#' Fits a shared parameter joint model for longitudinal and time-to-event 
#' (e.g. survival) data under a Bayesian framework using Stan.
#' 
#' @export
#' @template return-stanjm-object
#' @template see-also
#' @template args-jmpriors
#' @template args-prior_PD
#' @template args-adapt_delta
#' @template args-QR
#' 
#' @param formulaLong A two-sided linear formula object describing both the 
#'   fixed-effects and random-effects parts of the longitudinal submodel  
#'   (see \code{\link[lme4]{glmer}} for details). For a multivariate joint 
#'   model this should be a list of such formula objects, with each element
#'   of the list providing the formula for one of the longitudinal submodels.
#' @param dataLong A data frame containing the variables specified in
#'   \code{formulaLong}. If fitting a multivariate joint model, then this can
#'   be either a single data frame which contains the data/variables for all 
#'   the longitudinal submodels, or it can be a list of data frames where each
#'   element of the list provides the data for one of the longitudinal 
#'   submodels.
#' @param formulaEvent A two-sided formula object describing the time-to-event
#'   submodel. The left hand side of the formula should be a \code{Surv()} 
#'   object. See \code{\link[survival]{Surv}}.
#' @param dataEvent A data frame containing the variables specified in
#'   \code{formulaEvent}.
#' @param time_var A character string identifying the name of the variable 
#'   in \code{dataLong} which represents time.
#' @param family The family (and possibly also the link function) for the 
#'   longitudinal submodel. See \code{\link[lme4]{glmer}} for details.
#' @param assoc A character string or character vector specifying the joint
#'   model association structure. Possible association structures that can
#'   be used (described in the \strong{Details} section below) include: "null"; 
#'   "etavalue" (the default); "etaslope"; "muvalue"; or "shared_b". 
#'   By default, "shared_b" includes all random effects in the shared random 
#'   effects association structure, however, a subset of the random effects can 
#'   be chosen by specifying their indices between parentheses as a suffix, for 
#'   example, "shared_b(1)" or "shared_b(1:3)" or "shared_b(1,2,4)", and so on. 
#'   For a multivariate joint model, different association structures can be used for 
#'   each longitudinal submodel by specifying a list of character strings or character
#'   vectors, with each element of the list specifying the desired association 
#'   structure for one of the longitudinal submodels. Setting \code{assoc} equal to 
#'   \code{NULL} will fit a joint model with no association structure (equivalent  
#'   to fitting separate longitudinal and time-to-event models). See the
#'   \strong{Examples} section below for some examples.
#' @param base_haz A character string indicating which baseline hazard to use
#'   for the time-to-event submodel. Currently the only option allowed is 
#'   \code{"weibull"} (the default).
#' @param quadnodes A numeric scalar giving the number of nodes to use for 
#'   the Gauss-Kronrod quadrature. The quadrature is used to approximate the 
#'   integral over the cumulative hazard in the likelihood function. Options 
#'   are 7, 11 and 15 (the default).
#' @param subsetLong,subsetEvent Same as subset in \code{\link[stats]{glm}}.
#'   However, if fitting a multivariate joint model and a list of data frames 
#'   is provided in \code{dataLong} then a corresponding list of subsets 
#'   must be provided in \code{subsetLong}.
#' @param na.action,contrasts Same as \code{\link[stats]{glm}}, but rarely 
#'   specified.
#' @param weights,offset Same as \code{\link[stats]{glm}}. Not currently 
#'   allowed.
#' @param centreLong,centreEvent A logical specifying whether the predictor
#'   matrix for the longitudinal submodel(s) or event submodel should be 
#'   centred. 
#' @param init The method for generating the initial values for the MCMC.
#'   The default is "random" (see \code{\link[rstan]{stan}}), however, 
#'   \code{"model_based"} initial values are also available, which uses those  
#'   obtained from fitting separate longitudinal and time-to-event models prior  
#'   to fitting the joint model. Other possibilities for specifying \code{init}
#'   are those described for \code{\link[rstan]{stan}}.     
#' @param ... Further arguments passed to 
#'   \code{\link[rstan]{sampling}} (e.g. \code{iter}, \code{chains}, 
#'   \code{cores}, etc.).
#' @param prior_covariance Cannot be \code{NULL}; see \code{\link{decov}} for
#'   more information about the default arguments.
#'
#' @details The \code{stan_jm} function can be used to fit a joint model (also 
#'   known as a shared parameter model) for longitudinal and time-to-event data 
#'   under a Bayesian framework using Stan, but without needing to write the 
#'   model code or create the data list for Stan. 
#'   The joint model may be univariate (with only one longitudinal submodel) or
#'   multivariate (with more than one longitudinal submodel). \cr
#'   \cr 
#'   For the longitudinal submodel a generalised linear mixed model is assumed 
#'   with any of the \code{\link[stats]{family}} choices allowed by 
#'   \code{\link[lme4]{glmer}} (\strong{Note:} At this stage only gaussian 
#'   outcomes are currently implemented). For the event submodel a parametric
#'   proportional hazards model is assumed. Currently only a Weibull baseline 
#'   hazard is permitted. Time-varying covariates are allowed in both the 
#'   longitudinal and event submodels. 
#'   The association structure for the joint model can be based on any of the 
#'   following parameterisations: current value of the linear predictor in the 
#'   longitudinal submodel; current expected value in the longitudinal 
#'   submodel;  
#'   first derivative (slope) in the longitudinal submodel; first derivative 
#'   (slope) for the linear predictor in the longitudinal submodel; shared random 
#'   effects; no association structure (equivalent to fitting separate longitudinal 
#'   and event models). The association type is most easily specified using the 
#'   \code{\link{assoc}} function. \cr
#'   \cr
#'   Bayesian estimation is performed via MCMC. The Bayesian model includes 
#'   independent priors on the 
#'   regression coefficients for both the longitudinal and event submodels, 
#'   including the association parameter(s) (in much the same way as the
#'   regression parameters in \code{\link{stan_glm}}) and
#'   priors on the terms of a decomposition of the covariance matrices of the
#'   group-specific parameters (in the same way as \code{\link{stan_glmer}}). 
#'   See \code{\link{priors}} for more information about the priors distributions
#'   that are available. \cr
#'   \cr
#'   Gauss-Kronrod quadrature is used to numerically evaluate the integral  
#'   over the cumulative hazard in the likelihood function for the joint model.
#'   The accuracy of the numerical approximation can be controlled using the
#'   number of quadrature nodes, specified through the \code{quadnodes} 
#'   argument. Using a higher number of quadrature nodes will result in a more 
#'   accurate approximation.
#'    
#' @examples
#' #####
#' # Univariate joint model, with association structure based on the 
#' # current value of the linear predictor
#' f1 <- stan_jm(formulaLong = logBili ~ year + (1 | id), 
#'               dataLong = pbcLong,
#'               formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'               dataEvent = pbcSurv,
#'               time_var = "year",
#'               chains = 1, iter = 1000, warmup = 500)
#' summary(f1, digits = 3) 
#'         
#' #####
#' # Univariate joint model, with association structure based on the 
#' # current value of the linear predictor and shared random intercept
#' f2 <- stan_jm(formulaLong = logBili ~ year + (1 | id), 
#'               dataLong = pbcLong,
#'               formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'               dataEvent = pbcSurv,
#'               assoc = c("etavalue", "shared_b"),
#'               time_var = "year",
#'               chains = 1, iter = 1000, warmup = 500)
#' summary(f2, digits = 3)          
#' 
#' ######
#' # Multivariate joint model, with association structure based 
#' # on the current value of the linear predictor in each longitudinal 
#' # submodel and shared random intercept from the second longitudinal 
#' # submodel (which is the first random effect in that submodel
#' # and is therefore indexed the "-1" suffix in the code below)
#' mv1 <- stan_jm(formulaLong = list(
#'         logBili ~ year + (1 | id), 
#'         albumin ~ sex + year + (1 + year | id)),
#'         dataLong = pbcLong,
#'         formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'         dataEvent = pbcSurv,
#'         assoc = list("etavalue", c("etavalue", "shared_b(1)")), 
#'         time_var = "year", adapt_delta = 0.75,
#'         chains = 1, iter = 1000, warmup = 500, refresh = 25)
#' summary(mv1, digits = 3)
#' 
#' # To include both the random intercept and random slope in the shared 
#' # random effects association structure for the second longitudinal 
#' # submodel, we could specify the following
#' update(mv1, assoc = list("etavalue", c("etavalue", "shared_b"))
#' # which would be equivalent to    
#' update(mv1, assoc = list("etavalue", c("etavalue", "shared_b(1,2)"))                         
#' # or    
#' update(mv1, assoc = list("etavalue", c("etavalue", "shared_b(1:2)"))     
#' 
#' ######
#' # Multivariate joint model, estimated using multiple MCMC chains 
#' # run in parallel across all available PC cores
#' mv2 <- stan_jm(formulaLong = list(
#'         logBili ~ year + (1 | id), 
#'         albumin ~ sex + year + (1 +  year | id)),
#'         dataLong = pbcLong,
#'         formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'         dataEvent = pbcSurv,
#'         assoc = list("etavalue", c("etavalue", "shared_b(1)")),
#'         time_var = "year",
#'         chains = 3, iter = 1000, warmup = 500, refresh = 25,
#'         cores = parallel::detectCores())
#' summary(mv2, digits = 3)            
#'
#' @import data.table
#' @importFrom rstanarm normal student_t cauchy hs hs_plus decov
#' @importFrom survival Surv cluster
#' @importFrom Matrix Matrix t cBind bdiag
#' @importFrom JMbayes dns
stan_jm <- function(formulaLong, dataLong, 
                    formulaEvent, dataEvent, 
                    time_var, id_var, family = gaussian,
                    assoc = "etavalue",
                    base_haz = c("weibull", "piecewise", "splines"), 
                    df, quadnodes = 15, 
                    subsetLong, subsetEvent, 
                    na.action = getOption("na.action", "na.omit"),
                    weights, offset, contrasts,
                    centreLong = FALSE, centreEvent = FALSE, 
                    init = "random", ...,				          
					          priorLong = normal(), priorLong_intercept = normal(),
                    priorLong_ops = priorLong_options(),
                    priorEvent = normal(), priorEvent_intercept = normal(),
					          priorEvent_ops = priorEvent_options(),
					          priorAssoc = normal(),
					          priorAssoc_ops = priorAssoc_options(),
                    prior_covariance = decov(), prior_PD = FALSE,
                    adapt_delta = 0.75, QR = FALSE) {


  #=============================
  # Pre-processing of arguments
  #=============================  
  
  # Check for arguments not yet implemented
#  if (!missing(weights)) 
#    stop("Weights not yet supported by stan_jm")
  if (!missing(offset)) 
    stop("Offsets not yet supported by stan_jm")
#  algorithm <- match.arg(algorithm)
#  if (algorithm %in% c("meanfield", "fullrank"))
#    stop ("Meanfield and fullrank algorithms not yet implemented",
#          "for stan_jm")
  if (QR) 
    stop("QR decomposition not yet implemented for stan_jm")     
#  if ((init == "model_based") && any(unlist(c(centreLong, centreEvent)))) 
#    stop("Cannot use model based initial values when 'centreLong = TRUE'",
#         " or 'centreEvent = TRUE'.")  
  if (missing(id_var)) id_var <- NULL
  if (missing(subsetLong)) subsetLong <- NULL

  # Matched call
  call <- match.call(expand.dots = TRUE)    
  mc <- match.call(expand.dots = FALSE)
  mc$time_var <- mc$id_var <- 
    mc$assoc <- mc$base_haz <- mc$df <- mc$quadnodes <- 
    mc$centreLong <- mc$centreEvent <- mc$init <- NULL
  mc$priorLong <- mc$priorLong_intercept <- mc$priorLong_ops <- 
    mc$priorEvent <- mc$priorEvent_intercept <- mc$priorEvent_ops <-
    mc$priorAssoc <- mc$priorAssoc_ops <-
    mc$prior_covariance <- mc$prior_PD <- mc$algorithm <- mc$scale <- 
    mc$concentration <- mc$shape <-
    mc$adapt_delta <- mc$... <- mc$QR <- NULL
  
  # Create call for longitudinal submodel  
  y_mc <- mc
  y_mc <- strip_nms(y_mc, "Long") 
  y_mc$formulaEvent <- y_mc$dataEvent <- y_mc$subsetEvent <- NULL

  # Is formulaLong a list?
  if (is(eval(y_mc$formula), "formula")) { # number of long. markers
    formula_list <- FALSE
    M <- 1L
  } else if (is.list(eval(y_mc$formula))) {
    formula_list <- TRUE
    M <- length(eval(y_mc$formula))
    if (!all(sapply(eval(y_mc$formula), function(x)
      (is(x, "formula")))))
      stop("'formulaLong' should be a formula object or a list of formula ",
           "objects")               
  } else {
    stop("'formulaLong' should be a formula object or, for a multivariate ",
         "joint model, a list of formula objects with length equal to the ",
         "desired number of longitudinal markers")
  }
  if (M == 1L) cat("Univariate joint model specified")
  if (M > 1L)  cat("Multivariate joint model specified")
  
  # Is dataLong a list?
  if (is.data.frame(eval(y_mc$data))) {
    data_list <- FALSE
  } else if (is(eval(y_mc$data), "list")) {
    data_list <- TRUE
    if (length(eval(y_mc$data)) != M)
      stop("dataLong appears to be a list of the incorrect length")
  } else {
    stop("'dataLong' should be a data frame or possibly a list of data ",
         "frames. The latter is only required when fitting a multivariate ",
         "joint model using different data for each longitudinal submodel.")
  }
  
  # Is subset a list?
  if (is.null(y_mc$subset)) {
    y_subset_list <- NULL
  } else if (is(eval(y_mc$subset), "list")) {
    y_subset_list <- TRUE
    if (length(eval(y_mc$subset)) != M)
      stop("subsetLong appears to be a list of the incorrect length")
  } else if (is.vector(eval(y_mc$subset))) {
    y_subset_list <- FALSE
  } else {
    stop("'subsetLong' should be a vector or possibly a list of vectors. ",
         "The latter is only required if fitting a multivariate joint ",
         "model and using a different subset of data for each ",
         "longitudinal submodel.")
  }
  
  # Is weights a list?
  if (is.null(y_mc$weights)) {
    y_weights_list <- NULL
  } else if (is.list(eval(y_mc$weights))) {
    y_weights_list <- TRUE
    if (length(eval(y_mc$weights)) != M)
      stop("weights argument appears to be a list of the incorrect length")
  } else if (is.vector(eval(y_mc$weights))) {
    y_weights_list <- FALSE
  } else {
    stop("'weights' argument should be a numeric vector or possibly a list of ",
         "vectors. The latter is required if fitting a multivariate joint ",
         "model and using different weights for each longitudinal submodel. ",
         "If supplying a list then NULL should be specified for those submodels ",
         "that do not require weights.")
  }
  if (is.null(y_weights_list)) {
    y_weights <- rep(list(NULL), M)
  } else if (!y_weights_list) {
    y_weights <- list(weights)
    if (length(y_weights) < M) y_weights <- rep(y_weights, M)
  } else {
    y_weights <- weights
  }
  lapply(y_weights, rstanarm:::validate_weights)
  
  # Is family a list?
  if (is.null(eval(y_mc$family))) {
    family_list <- NULL
  } else if (is(eval(y_mc$family), "list")) {
    family_list <- TRUE
    if (length(eval(y_mc$family)) != M)
      stop("'family' argument must be a family or possibly a list of families. ",
           "The latter is only required when fitting a multivariate joint ",
           "model with a different family and/or link function for some ",
           "of the longitudinal submodels.")
  } else family_list <- FALSE
  if (is.null(family_list) || !family_list) family <- list(family)  # convert to list
  if (length(family) < M) family <- rep(family, M)  # repeat family if necessary
  
  # Is assoc a list? If not, then convert to list
  if (is.list(assoc)) {  # if list, then check length
    if (!(length(assoc) %in% c(1,M)))
      stop("`assoc' should be a list of length 1 or length equal to the ",
           "number of longitudinal markers")
  } else if (is.character(assoc)) {  # if not list, then convert to list
    assoc <- list(assoc)
  } else {  # else return error
    stop("'assoc' should be a character string or character vector or, for a ",
         "multivariate joint model, possibly a list of character strings ",
         "or character vectors. The latter is only required if using a different ",
         "association structure for linking each longitudinal submodel to the ",
         "event outcome.")
  }
  if (length(assoc) != M) assoc <- rep(assoc, M)
  
  # Check family and link
  supported_families <- c("binomial", "gaussian", "Gamma", "inverse.gaussian",
                          "poisson", "neg_binomial_2")
  family <- lapply(family, rstanarm:::validate_family)
  fam <- lapply(family, function(x) 
                which(pmatch(supported_families, x$family, nomatch = 0L) == 1L))
  if (any(lapply(fam, length) == 0L)) 
    stop("'family' must be one of ", paste(supported_families, collapse = ", "))
  supported_links <- lapply(fam, function(x) 
    switch(
      supported_families[x],
      binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
      gaussian = c("identity", "log", "inverse"),
      Gamma = c("identity", "log", "inverse"),
      inverse.gaussian = c("identity", "log", "inverse", "1/mu^2"),
      "neg_binomial_2" = , # intentional
      poisson = c("log", "identity", "sqrt"),
      stop("unsupported family")
    )
  )
  link <- mapply(function(x, i) which(supported_links[[i]] == x$link),
                 family, seq_along(family), SIMPLIFY = TRUE)
  if (any(lapply(link, length) == 0L)) 
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))

  # Create call for each longitudinal submodel separately
  m_mc <- list()  # list containing matched calls for each marker
  for (m in 1:M) {
    m_mc[[m]] <- y_mc
    m_mc[[m]]$formula <- if (formula_list) y_mc$formula[[(1+m)]] else y_mc$formula
    m_mc[[m]]$data    <- if (data_list)    y_mc$data[[(1+m)]]    else y_mc$data
    if (!is.null(y_subset_list))   
      m_mc[[m]]$subset  <- if (y_subset_list) y_mc$subset[[(1+m)]]  else y_mc$subset
    if (!is.null(y_weights_list))   
      m_mc[[m]]$weights  <- if (y_weights_list) y_mc$weights[[(1+m)]]  else y_mc$weights
    if (!is.null(family_list))   
      m_mc[[m]]$family  <- if (family_list)   y_mc$family[[(1+m)]]  else y_mc$family    
  }
   
  # Create call for event submodel
  e_mc <- mc
  e_mc <- strip_nms(e_mc, "Event")
  e_mc$formulaLong <- e_mc$dataLong <- e_mc$family <-
    e_mc$subsetLong <- NULL
                          
  # Standardised GK quadrature points and weights
  quadpoints <- get_quadpoints(quadnodes)
            
     
  #================================
  # Data for longitudinal submodel
  #================================
  
  # Items to store for each longitudinal submodel
  y_mod         <- list()     # fitted long. submodels
  y_vars        <- list()     # the variables used in fitting the model
  y_is_real     <- c()        # indicator of response vector being reals (not integers)
  y             <- list()     # response vector
  x             <- list()     # design matrix with intercept
  xtemp         <- list()     # design matrix without intercept, possibly centred
  xbar          <- list()     # means of predictors
  trials        <- list()     # num. of trials for binomial outcome
  y_centre      <- c()        # submodel has intercept
  y_has_intercept <- c()      # submodel has intercept
  y_has_intercept_unbound <- c()    # has unbounded intercept
  y_has_intercept_lobound <- c()    # has lower bounded intercept
  y_has_intercept_upbound <- c()    # has upper bounded intercept
  y_has_weights <- c()        # submodel has weights
  y_has_dispersion <- c()     # submodel has dispersion term
  y_N           <- c()        # num. observations
  y_N01         <- list()     # num. 0 and 1 observations if bernoulli
  y_real_N      <- c()        # num. observations, for real outcomes
  y_int_N       <- c()        # num. observations, for integer outcomes
  y_K           <- c()        # num. predictors (excluding intercept)
  y_offset      <- list()     # offsets
  Z             <- list()     # Z matrices
  y_cnms          <- list()   
  y_flist         <- list()   
  y_gamma_unbound <- list()   # initial values for intercepts
  y_gamma_lobound <- list()   # initial values for intercepts
  y_gamma_upbound <- list()   # initial values for intercepts
  y_beta          <- list()   # initial values for coefs
  y_dispersion    <- c()      # initial values for dispersion
  sd_b            <- list()   # initial values for random effect sds
  b_Corr          <- list()   # initial values for correlation matrix

  for (m in 1:M) {
  
    #if (M == 1) 
      #cat("\n--> Fitting separate longitudinal model now...") 
    #else 
      #cat(paste0("\n--> Fitting separate longitudinal model for marker ", m, " now..."))  
      
    # Fit separate longitudinal model
    if ((family[[m]]$family == "gaussian") && (family[[m]]$link == "identity")) {
      m_mc[[m]][[1]] <- quote(lme4::lmer)
      m_mc[[m]]$family <- NULL
      m_mc[[m]]$control <- get_control_args()
    } else {
      m_mc[[m]][[1]] <- quote(lme4::glmer)
      m_mc[[m]]$control <- get_control_args(glmer = TRUE)               
    }
    y_mod[[m]] <- eval(m_mc[[m]], parent.frame())      

    # Indicator of real or integer response vector
    y_is_real[m] <- check_response_real(family[[m]]$family)
    
    # Indicator of whether model includes a dispersion term
    y_has_dispersion[m] <- check_for_dispersion(family[[m]]$family)
    
    # Response vector and design matrix
    y[[m]] <- as.vector(lme4::getME(y_mod[[m]], "y"))
    x[[m]] <- as.matrix(lme4::getME(y_mod[[m]], "X"))
    y_has_intercept[m] <- grepl("(Intercept", colnames(x[[m]])[1L], fixed = TRUE)
    if (y_has_intercept[m]) {
      check_int <- check_intercept(family[[m]]$family, family[[m]]$link)
      y_has_intercept_unbound[m] <- check_int$unbound
      y_has_intercept_lobound[m] <- check_int$lobound
      y_has_intercept_upbound[m] <- check_int$upbound
    } else {
      y_has_intercept_unbound[m] <- 0L
      y_has_intercept_lobound[m] <- 0L
      y_has_intercept_upbound[m] <- 0L
    } 
    xtemp[[m]] <- if (y_has_intercept[m]) x[[m]][, -1L, drop=FALSE] else x[[m]]

    if (rstanarm:::is.binomial(family[[m]]$family)) {
      if (NCOL(y[[m]]) == 1L) {
        if (is.numeric(y[[m]]) || is.logical(y[[m]])) 
          y[[m]] <- as.integer(y[[m]])
        if (is.factor(y[[m]])) 
          y[[m]] <- fac2bin(y[[m]])
        if (!all(y[[m]] %in% c(0L, 1L))) 
          stop("y values must be 0 or 1 for bernoulli regression.")
        trials[[m]] <- rep(0L, length(y[[m]]))
      } else {
        if (!isTRUE(NCOL(y[[m]]) == 2L))
          stop("y should either be a vector or a matrix 1 or 2 columns.")
        trials[[m]] <- as.integer(y[[m]][, 1L] + y[[m]][, 2L])
        y[[m]] <- as.integer(y[[m]][, 1L])
      }
    } else trials[[m]] <- rep(0L, length(y[[m]]))

    # Falsify weights vector if weights are not included
    if (!length(y_weights[[m]])) {
      y_weights[[m]] <- rep(0.0, length(y[[m]]))
      y_has_weights[m] <- 0L
    } else y_has_weights[m] <- 1L
    
    # Random effect terms
    Z[[m]]     <- lme4::getME(y_mod[[m]], "Z")
    y_cnms[[m]]  <- lme4::getME(y_mod[[m]], "cnms")
    y_flist[[m]] <- lme4::getME(y_mod[[m]], "flist")
    y_offset[[m]] <- lme4::getME(y_mod[[m]], "offset")
        
    # Centred design matrix, if required
    if (centreLong) {
      y_centre[m] <- 1L
      xbar[[m]] <- colMeans(xtemp[[m]])
      xtemp[[m]] <- sweep(xtemp[[m]], 2, xbar[[m]], FUN = "-")
    }
    
    # Reorder y, X, Z if bernoulli (zeros first)
    if (rstanarm:::is.binomial(family[[m]]$family) && all(y[[m]] %in% 0:1)) {      
      y[[m]] <- y[[m]][order(y[[m]])]
      trials[[m]] <- trials[[m]][order(y[[m]])]
      y_weights[[m]] <- y_weights[[m]][order(y[[m]])]
      xtemp[[m]] <- xtemp[[m]][order(y[[m]]), , drop = FALSE]  
      Z[[m]] <- Z[[m]][order(y[[m]]), , drop = FALSE]
      y_N01[[m]] <- sapply(0:1, function(x) sum(y[[m]] == x))
    } else y_N01[[m]] <- rep(0L, 2)  # dud entry if not bernoulli
    
    # Dimensions
    y_N[m] <- NROW(xtemp[[m]])
    y_real_N[m] <- if (y_is_real[m]) y_N[m] else 0L
    y_int_N[m] <- if (!y_is_real[m]) y_N[m] else 0L
    y_K[m] <- NCOL(xtemp[[m]])

    # Update formula if using splines or other data dependent predictors
    y_vars[[m]] <- get_formvars(y_mod[[m]])

    # Model based initial values
    if (init == "model_based") {
      if (y_has_intercept_unbound[m]) {
        y_beta[[m]] <- lme4::fixef(y_mod[[m]])[-1L]
        markun <- sum(y_has_intercept_unbound[1:m])
        y_gamma_unbound[[markun]] <- fixef(y_mod[[m]])[1L]
        if (centreLong) 
          y_gamma_unbound[[markun]] <- y_gamma_unbound[[markun]] - xbar[[m]] %*% y_beta[[m]]
      } else if (y_has_intercept_lobound[m]) {
        y_beta[[m]] <- lme4::fixef(y_mod[[m]])[-1L]
        marklo <- sum(y_has_intercept_lobound[1:m])        
        y_gamma_lobound[[marklo]] <- fixef(y_mod[[m]])[1L]
        if (centreLong) 
          y_gamma_lobound[[marklo]] <- y_gamma_lobound[[marklo]] - xbar[[m]] %*% y_beta[[m]]
      } else if (y_has_intercept_upbound[m]) {
        y_beta[[m]] <- lme4::fixef(y_mod[[m]])[-1L]
        markup <- sum(y_has_intercept_upbound[1:m])        
        y_gamma_upbound[[markup]] <- fixef(y_mod[[m]])[1L]
        if (centreLong) 
          y_gamma_upbound[[markup]] <- y_gamma_upbound[[markup]] - xbar[[m]] %*% y_beta[[m]]
      } else {
        y_beta[[m]] <- lme4::fixef(y_mod[[m]])
      }
      vc <- lme4::VarCorr(y_mod[[m]])[[1]]
      sd_b[[m]] <- attr(vc, "stddev")
      b_Corr[[m]] <- attr(vc, "correlation")
      if (y_has_dispersion[m]) {
        disp_mark <- sum(y_has_dispersion[1:m])
        y_dispersion[disp_mark] <- sigma(y_mod[[m]])
      }
    }

  }
  
  # Sum dimensions across all longitudinal submodels
  sum_y_N <- sum(y_N)
  sum_y_real_N <- sum(y_real_N)
  sum_y_int_N <- sum(y_int_N)
  sum_y_K <- sum(y_K)
  sum_y_has_intercept <- sum(y_has_intercept)
  sum_y_has_intercept_unbound <- sum(y_has_intercept_unbound)
  sum_y_has_intercept_lobound <- sum(y_has_intercept_lobound)
  sum_y_has_intercept_upbound <- sum(y_has_intercept_upbound)
  sum_y_has_dispersion <- sum(y_has_dispersion)
  
  # Indexing for binded response vector, design matrix, weights, etc
  y_beg <- sapply(1:M, function(m) sum(y_N[0:(m-1)]) + 1)
  y_end <- sapply(1:M, function(m) sum(y_N[0:m]))  
  y_real_beg <- sapply(1:M, function(m) if (y_is_real[m]) sum(y_real_N[0:(m-1)]) + 1L else 0L)
  y_real_end <- sapply(1:M, function(m) if (y_is_real[m]) sum(y_real_N[0:m]) else 0L)  
  y_int_beg <- sapply(1:M, function(m) if (!y_is_real[m]) sum(y_int_N[0:(m-1)]) + 1L else 0L)
  y_int_end <- sapply(1:M, function(m) if (!y_is_real[m]) sum(y_int_N[0:m]) else 0L)  
  
  # Additional error checks
  id_var <- check_id_var(id_var, y_cnms)
  id_list <- check_id_list(id_var, y_flist)
  
  # Construct single cnms list for all longitudinal submodels
  y_cnms_nms <- lapply(y_cnms, names)
  cnms_nms <- unique(unlist(y_cnms_nms))
  cnms <- lapply(seq_along(cnms_nms), function(i) {
    nm <- cnms_nms[i]
    unlist(lapply(1:length(y_cnms), function(m) 
      if (nm %in% y_cnms_nms[[m]]) paste0("Long", m, "|", y_cnms[[m]][[nm]])))
  })
  names(cnms) <- cnms_nms
  
  # Family indicators
  famname <- lapply(fam, function(x) supported_families[x])
  is_bernoulli  <- mapply(function(x, i)
    rstanarm:::is.binomial(x) && all(y[[i]] %in% 0:1),
                          famname, seq_along(famname), SIMPLIFY = FALSE)
  is_nb         <- lapply(famname, rstanarm:::is.nb)
  is_gaussian   <- lapply(famname, rstanarm:::is.gaussian)
  is_gamma      <- lapply(famname, rstanarm:::is.gamma)
  is_ig         <- lapply(famname, rstanarm:::is.ig)
  is_continuous <- lapply(seq_along(famname), function(x) 
                     (is_gaussian[[x]] || is_gamma[[x]] || is_ig[[x]]))
  
  # Require intercept for certain family and link combinations
  lapply(1:M, function(x) {
    if (!y_has_intercept[x]) {
      linkname <- supported_links[[x]][link[[x]]]
      needs_intercept <- 
        !is_gaussian[[x]] && linkname == "identity" ||
        is_gamma[[x]] && linkname == "inverse" ||
        rstanarm:::is.binomial(famname[[x]]) && linkname == "log"
      if (needs_intercept)
        stop(paste0("To use the combination of family and link ", 
                    "specified for longitudinal marker ", x,
                    ", the model must have an intercept."))
    }
  })
    
  
  #=========================
  # Data for event submodel
  #=========================

  # Items to store for event submodel
  e_beta <- c()
  
  # Survival submodel
  #cat("\n--> Fitting separate survival model now...") 
  
  # Baseline hazard
  base_haz <- match.arg(base_haz)
  base_haz_weibull <- (base_haz == "weibull")
  base_haz_piecewise <- (base_haz == "piecewise")
  base_haz_splines <- (base_haz == "splines")
  
  if ((base_haz_splines) && (missing(df)))
    stop("df must be specified when splines baseline hazard is used",
         call. = FALSE)
  if ((!base_haz_splines) && (!missing(df)))
    stop("df can only be specified when splines baseline hazard is used",
         call. = FALSE)
  if (missing(df)) {
    df <- NULL
    splines_df <- 0L
  } else splines_df <- df
  
  # Set up model frame for event submodel 
  e_mc[[1]] <- quote(survival::coxph) 
  e_mc$x <- TRUE
  e_mod <- eval(e_mc, parent.frame())
  e_fr <- e_mf <- expand.model.frame(e_mod, id_var, na.expand = TRUE)
  e_mf <- cbind(unclass(e_mf[,1]), e_mf)

  # Check ID sorting
  e_id_list <- factor(unique(e_mf[, id_var]))
  if (!identical(id_list, e_id_list))
    stop("'dataEvent' needs to be sorted by the subject ",
         "ID/grouping variable", call. = FALSE)

  e_y <- e_mod$y
  
  # For each individual, identify final event time and event indicator
  if (attr(e_y, "type") == "counting") {
    tvc         <- TRUE
    mf_event    <- do.call(rbind, lapply(
                             split(e_mf, e_mf[, id_var]),
                             function(d) d[which.max(d[,"stop"]), ]))
    flist_event <- mf_event[, id_var]
    eventtime   <- mf_event$stop
    d           <- mf_event$status
  
    e_mf           <- data.table(cbind(e_y, e_mf), key = c(id_var, "start"))
    e_mf_eventtime <- e_mf[, .SD[.N], by = e_mf[, id_var]]
    # Unstandardised quadrature points
    quadpoint <- lapply(quadpoints$points, FUN = function(x) 
                          (eventtime/2) * x + (eventtime/2))
    # Model frame corresponding to observation times which are 
    #   as close as possible to the unstandardised quadrature points                      
    e_mf_quadtime  <- do.call(rbind, lapply(quadpoint, FUN = function(x) 
                         e_mf[data.table::SJ(flist_event, x), 
                         roll = TRUE, rollends = c(TRUE, TRUE)]))
    # Model frame evaluated at both event times and quadrature points
    e_mf_quadtime <- rbind(e_mf_eventtime, e_mf_quadtime, idcol = "xbind.id")
    # Design matrix evaluated at event times and quadrature points
    #   NB Here there are time varying covariates in the event submodel and
    #   therefore the design matrix differs depending on the quadrature point 
    e_x_quadtime   <- update(e_mod, data = e_mf_quadtime)$x
  } else if (attr(e_y, "type") == "right") {
    tvc         <- FALSE 
    mf_event    <- e_mf
    flist_event <- mf_event[, id_var]
    eventtime   <- mf_event$time
    d           <- mf_event$status
    # Unstandardised quadrature points
    quadpoint <- lapply(quadpoints$points, FUN = function(x) 
                          (eventtime/2) * x + (eventtime/2))    
    # Design matrix evaluated at event times and quadrature points
    #   NB Here there are no time varying covariates in the event submodel and
    #   therefore the design matrix is identical at event time and at all
    #   quadrature points
    e_x_quadtime   <- do.call(rbind, lapply(1:(quadnodes + 1), 
                                            FUN = function(x) e_mod$x))
  } else stop("Only 'right' or 'counting' type Surv objects are allowed 
               on the LHS of the event submodel formula")

  # Incorporate intercept term (since Cox model does not have intercept)
  e_x_quadtime  <- cbind("(Intercept)" = rep(1, NROW(e_x_quadtime)), e_x_quadtime)

  # Centering of design matrix for event model
  e_x <- as.matrix(e_x_quadtime)  
  e_has_intercept <- grepl("(Intercept", colnames(e_x)[1L], fixed = TRUE)
  e_xtemp <- if (e_has_intercept) e_x[, -1L, drop=FALSE] else e_x
  if (centreEvent) {
    e_xbar <- colMeans(e_xtemp)
    e_xtemp <- sweep(e_xtemp, 2, e_xbar, FUN = "-")
  }
  e_K <- NCOL(e_xtemp)
  Npat <- length(eventtime)
  weights_rep <- rep(quadpoints$weights, each = Npat)  
  eventtime_rep <- rep(eventtime, times = quadnodes)  
  weights_times_half_eventtime <- 0.5 * weights_rep * eventtime_rep   
  
  # Model based initial values
  if (init == "model_based") {  
    e_beta <- e_mod$coef
  }
    
  # Error checks for the ID variable
  id_list <<- id_list
  flist_event <<- flist_event
  if (!identical(id_list, factor(sort(unique(flist_event)))))
    stop("The patient IDs (levels of the grouping factor) included ",
         "in the longitudinal and event submodels do not match")


  #================================
  # Data for association structure
  #================================
  
  # Check association structure
  supported_assocs <- c("null", "etavalue", "etaslope", "muvalue", "muslope", "shared_b")
  assoc <- lapply(assoc, validate_assoc, supported_assocs)

  # Indicator of each association type, for each longitudinal submodel
  has_assoc <- sapply(supported_assocs, function(x) 
    sapply(assoc, function(y) as.integer(y[[x]])), simplify = FALSE)
  
  # Shared random effects
  which_b <- lapply(1:M, function(m) assoc[[m]]$which_b)
  if (any(has_assoc$shared_b)) {
    max_which_b <- sapply(y_cnms, function(x) length(x[[id_var]]))
    for (m in 1:M) {
      if ((has_assoc$shared_b[m]) && (!length(which_b[[m]]))) 
        which_b[[m]] <- seq_len(max_which_b[m])
      if (any(which_b[[m]] > max_which_b[m]))
        stop(paste0("The indices specified for the shared random effects association ",
                    "structure are greater than the number of subject-specific random ", 
                    "effects (this error was encountered for longitudinal submodel ", m,
                    ")"), call. = FALSE)
    }
  }
  size_which_b <- sapply(which_b, length)
  a_K <- get_num_assoc_pars(has_assoc, which_b)


  #====================================================================
  # Longitudinal submodel: calculate design matrices, and id vector 
  # at the event times and quadrature points
  #====================================================================
    
  # Items to store for each longitudinal submodel
  y_mod_q         <- list()   # fitted long. submodels at quadpoints
  y_fr            <- list()   # model frame with the addition of time_var
  xq              <- list()   # design matrix before removing intercept 
                              # and centering
  xqtemp          <- list()   # design matrix (without intercept) for 
                              # longitudinal submodel calculated at event 
                              # and quad times, possibly centred
  dxdtq           <- list()   # first derivative of design matrix
  dxdtqtemp       <- list()   # first derivative of design matrix (intercept 
                              # term removed)
  Zq              <- list()   # random effects matrix
  dZdtq           <- list()   # first derivative of random effects matrix

  # Set up a second longitudinal model frame which includes the time variable
  for (m in 1:M) {
    m_mc_temp         <- m_mc[[m]]
    m_mc_temp[[1]]    <- quote(lme4::glFormula)    
    m_mc_temp$formula <- do.call(update.formula, list(
                                    m_mc[[m]]$formula, 
                                    paste0("~ . +", time_var))) 
    if ((family[[m]]$family == "gaussian") && (family[[m]]$link == "identity")) {
      m_mc_temp$control <- get_control_args(norank = TRUE)
    } else {
      m_mc_temp$control <- get_control_args(glmer = TRUE, norank = TRUE)
    }     
    
    # Obtain model frame with time_var definitely included
    y_mod_wtime <- eval(m_mc_temp, parent.frame())      
    y_fr[[m]] <- y_mod_wtime$fr
    
    # Update model based on predvars and obtain new model frame
    mf <- data.table(y_mod_wtime$fr, key = c(id_var, time_var))
  
    # Identify which row in longitudinal data is closest to event time
    mf_eventtime <- mf[data.table::SJ(flist_event, eventtime), 
                          roll = TRUE, rollends = c(FALSE, TRUE)]
  
    # Identify which row in longitudinal data is closest to quadrature point
    #   NB if the quadrature point is earlier than the first observation time, 
    #   then covariates values are carried back to avoid missing values - I
    #   should add a warning for when this is the case! In any other case, the 
    #   observed covariates values from the most recent observation time
    #   preceeding the quadrature point are carried forward to represent the 
    #   covariate value(s) at the quadrature point. (To avoid missingness  
    #   there is no limit on how far forwards or how far backwards covariate 
    #   values can be carried). If no time varying covariates are present in
    #   the longitudinal submodel (other than the time variable) then nothing 
    #   is carried forward or backward.
    mf_quadtime <- do.call(rbind, lapply(quadpoint, FUN = function(x) 
                           mf[data.table::SJ(flist_event, x), 
                           roll = TRUE, rollends = c(TRUE, TRUE)]))
  
    # Obtain long design matrix evaluated at event times (xbind.id == 1) and  
    #   quadrature points (xbind.id == 2)
    names(mf_eventtime)[names(mf_eventtime) == "eventtime"] <- time_var
    names(mf_quadtime)[names(mf_quadtime) == "quadpoint"]   <- time_var
    mf_quadtime <- rbind(mf_eventtime, mf_quadtime, idcol = "xbind.id")
    m_mc_temp$formula <- m_mc[[m]]$formula  # return to original formula
    m_mc_temp$formula <- substitute_vars(y_mod[[m]], y_vars[[m]]$formvars, y_vars[[m]]$predvars)
    m_mc_temp$control <- m_mc[[m]]$control  # return to original control args
    m_mc_temp$data    <- mf_quadtime        # data at event and quadrature times
    y_mod_q[[m]] <- eval(m_mc_temp, parent.frame())     
         
    xq[[m]] <- as.matrix(y_mod_q[[m]]$X)
    xqtemp[[m]] <- if (y_has_intercept[m]) xq[[m]][, -1L, drop=FALSE] else xq[[m]]  
    Zq[[m]] <- t(y_mod_q[[m]]$reTrms$Zt)
      #Needs working out to appropriately deal with offsets??
      #offset_quadtime <- model.offset(mod_quadtime$fr) %ORifNULL% double(0)

    # Centering of design matrix for longitudinal model at event times
    # and quadrature times 
    if (centreLong) xqtemp[[m]] <- sweep(xqtemp[[m]], 2, xbar[[m]], FUN = "-")
  }

  if(sum(has_assoc$etaslope))
    cat("\nSlope association structure based on the following first derivatives:")
  for (m in 1:M) {
    if (has_assoc$etaslope[m]) {
      sel_lin <- lapply(y_vars[[m]]$predvars, function(x) 
                       grep(paste0("^", time_var, "$"), x, value = FALSE))
      sel_ns <- lapply(y_vars[[m]]$predvars, function(x) 
                       grep(paste0("^ns\\(.*", time_var, ".*\\)"), x, value = FALSE))
      # Replace those variables with their first derivative
      formvars_deriv <- y_vars[[m]]$formvars      
      predvars_deriv <- y_vars[[m]]$predvars
      for (i in 1:length(sel_lin)) {
        formvars_deriv$fixed[sel_lin$fixed] <- 
          gsub(paste0("^", time_var, "$"), "1", formvars_deriv$fixed[sel_lin$fixed])
        formvars_deriv$random[sel_lin$random] <- 
          gsub(paste0("^", time_var, "$"), "1", formvars_deriv$random[sel_lin$random])
        predvars_deriv$fixed[sel_lin$fixed] <- 
          gsub(paste0("^", time_var, "$"), "1", predvars_deriv$fixed[sel_lin$fixed])
        predvars_deriv$random[sel_lin$random] <- 
          gsub(paste0("^", time_var, "$"), "1", predvars_deriv$random[sel_lin$random])
      }
      for (i in 1:length(sel_ns)) {
        formvars_deriv$fixed[sel_ns$fixed] <- 
          gsub("^ns\\(", "dns(", formvars_deriv$fixed[sel_ns$fixed])
        formvars_deriv$random[sel_ns$random] <- 
          gsub("^ns\\(", "dns(", formvars_deriv$random[sel_ns$random])
        predvars_deriv$fixed[sel_ns$fixed] <- 
          gsub("^ns\\(", "dns(", predvars_deriv$fixed[sel_ns$fixed])
        predvars_deriv$random[sel_ns$random] <- 
          gsub("^ns\\(", "dns(", predvars_deriv$random[sel_ns$random])      
      }      
      # Tell user which variables have been replaced
      cat(paste0("\n  Submodel ", m, ": "), 
          paste0(y_vars[[m]]$formvars$fixed[c(sel_lin$fixed, sel_ns$fixed)], " -> ", 
                 formvars_deriv$fixed[c(sel_lin$fixed, sel_ns$fixed)], collapse = "; "),
          paste0(y_vars[[m]]$formvars$random[c(sel_lin$random, sel_ns$random)], " -> ", 
                 formvars_deriv$random[c(sel_lin$random, sel_ns$random)], collapse = "; "))
      # Replace variables in model formula
      m_mc_temp$formula <- substitute_vars(y_mod_q[[m]], y_vars[[m]]$predvars, predvars_deriv)
      if ((family[[m]]$family == "gaussian") && (family[[m]]$link == "identity")) {
        m_mc_temp$control <- get_control_args(norank = TRUE)
      } else {
        m_mc_temp$control <- get_control_args(glmer = TRUE, norank = TRUE)      
      }
      # Evaluate X and Z using derivatives
      y_mod_q_dxdt <- eval(m_mc_temp, parent.frame())
      # Obtain new X matrix
      dxdtq[[m]] <- as.matrix(y_mod_q_dxdt$X)
      dxdtqtemp[[m]] <- if (y_has_intercept[m]) dxdtq[[m]][, -1L, drop=FALSE] else dxdtq[[m]]
      # Identify columns in X which have a non-zero derivative and replace all
      # other columns with a value of zero
      sel <- grepl(paste0("dns\\(.*", time_var, ".*\\)"), colnames(dxdtqtemp[[m]]))
      dxdtqtemp[[m]][, !sel] <- 0
      # Obtain new Z matrix
      dZdtq[[m]] <- t(y_mod_q_dxdt$reTrms$Zt)
      # For each grouping factor, identify which terms involve a non-zero derivative
      sel <- lapply(y_mod_q_dxdt$reTrms$cnms, function(x) grepl(paste0("dns\\(.*", time_var, ".*\\)"), x))
      # Calculate number of factor levels for each grouping factor
      len_flist <- sapply(y_mod_q_dxdt$reTrms$flist, function(x) length(unique(x)))
      # Using the selection of terms, and the number of factor levels, create a 
      # vector indicating which columns of Z have a non-zero derivative and 
      # replace all other columns with a value of zero
      sel_vec <- unlist(lapply(seq_along(sel), function(x) rep(sel[[x]], len_flist[x])))
      dZdtq[[m]][, !sel_vec] <- 0 
    } else {
      dxdtqtemp[[m]] <- matrix(0, ncol = ncol(xqtemp[[m]]), nrow = nrow(xqtemp[[m]]))
      dZdtq[[m]] <- Matrix::sparseMatrix(dims = c(nrow(Zq[[m]]),ncol(Zq[[m]])), i={}, j={}, x = 0)
    }
  }
  if(sum(has_assoc$etaslope))
    cat("\nAny interactions involving these variables will also be preserved.")


  #=====================
  # Prior distributions
  #=====================
 
  ok_dists <- rstanarm:::nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus")
  ok_intercept_dists <- ok_dists[1:3]

  # Priors for longitudinal submodel(s)
  priorLong_scaled <- priorLong_ops$scaled
  priorLong_min_prior_scale <- priorLong_ops$min_prior_scale
  priorLong_scale_for_dispersion <- 
    as.array(rstanarm:::maybe_broadcast(priorLong_ops$prior_scale_for_dispersion, sum_y_has_dispersion))
  
  if (is.null(priorLong)) {
    priorLong_dist <- 0L
    priorLong_mean <- as.array(rep(0, sum_y_K))
    priorLong_scale <- priorLong_df <- as.array(rep(1, sum_y_K))
  } else {
    if (!is.list(priorLong)) 
      stop("'prior' should be a named list.")
    priorLong_dist <- priorLong$dist
    priorLong_scale <- priorLong$scale
    priorLong_mean <- priorLong$location
    priorLong_df <- priorLong$df
    priorLong_df[is.na(priorLong_df)] <- 1
    if (!priorLong_dist %in% unlist(ok_dists)) {
      stop("The prior distribution for the coefficients should be one of ",
           paste(names(ok_dists), collapse = ", "))
    } else if (priorLong_dist %in% c("normal", "t")) {
      priorLong_dist <- ifelse(priorLong_dist == "normal", 1L, 2L)
      # !!! Need to change this so that link can be specific to each submodel
      priorLong_scale <- rstanarm:::set_prior_scale(priorLong_scale, default = 2.5, 
                                     link = family[[1]]$link)
    } else {
      priorLong_dist <- ifelse(priorLong_dist == "hs", 3L, 4L)
    }
    
    priorLong_df <- rstanarm:::maybe_broadcast(priorLong_df, sum_y_K)
    priorLong_df <- as.array(pmin(.Machine$double.xmax, priorLong_df))
    priorLong_mean <- rstanarm:::maybe_broadcast(priorLong_mean, sum_y_K)
    priorLong_mean <- as.array(priorLong_mean)
    priorLong_scale <- rstanarm:::maybe_broadcast(priorLong_scale, sum_y_K)
    priorLong_scale <- as.array(priorLong_scale)
  }
  if (is.null(priorLong_intercept)) {
    priorLong_dist_for_intercept <- 0L
    priorLong_mean_for_intercept <- as.array(rep(0, sum_y_has_intercept))
    priorLong_scale_for_intercept <- priorLong_df_for_intercept <- as.array(rep(1, sum_y_has_intercept))
  } else {
    if (!is.list(priorLong_intercept)) 
      stop("'priorLong_intercept' should be a named list.")
    priorLong_dist_for_intercept <- priorLong_intercept$dist
    priorLong_scale_for_intercept <- priorLong_intercept$scale
    priorLong_mean_for_intercept <- priorLong_intercept$location
    priorLong_df_for_intercept <- priorLong_intercept$df 
    priorLong_df_for_intercept[is.na(priorLong_df_for_intercept)] <- 1
    
    if (!priorLong_dist_for_intercept %in% unlist(ok_intercept_dists))
      stop("The prior distribution for the intercept should be one of ",
           paste(names(ok_intercept_dists), collapse = ", "))
    priorLong_dist_for_intercept <- 
      ifelse(priorLong_dist_for_intercept == "normal", 1L, 2L)
    priorLong_scale_for_intercept <- 
      rstanarm:::set_prior_scale(priorLong_scale_for_intercept, default = 10, 
                      link = family[[1]]$link)

    priorLong_df_for_intercept <- rstanarm:::maybe_broadcast(priorLong_df_for_intercept, sum_y_has_intercept)
    priorLong_df_for_intercept <- as.array(pmin(.Machine$double.xmax, priorLong_df_for_intercept))
    priorLong_mean_for_intercept <- rstanarm:::maybe_broadcast(priorLong_mean_for_intercept, sum_y_has_intercept)
    priorLong_mean_for_intercept <- as.array(priorLong_mean_for_intercept)
    priorLong_scale_for_intercept <- rstanarm:::maybe_broadcast(priorLong_scale_for_intercept, sum_y_has_intercept)
    priorLong_scale_for_intercept <- as.array(priorLong_scale_for_intercept)
  }

  # Priors for event submodel
  priorEvent_scaled <- priorEvent_ops$scaled
  priorEvent_min_prior_scale <- priorEvent_ops$min_prior_scale
  priorEvent_scale_for_weibull <- priorEvent_ops$prior_scale_for_weibull
  
  if (is.null(priorEvent)) {
    priorEvent_dist <- 0L
    priorEvent_mean <- as.array(rep(0, e_K))
    priorEvent_scale <- priorEvent_df <- as.array(rep(1, e_K))
  } else {
    if (!is.list(priorEvent)) 
      stop("'priorEvent' should be a named list.")
    priorEvent_dist <- priorEvent$dist
    priorEvent_scale <- priorEvent$scale
    priorEvent_mean <- priorEvent$location
    priorEvent_df <- priorEvent$df
    priorEvent_df[is.na(priorEvent_df)] <- 1
    if (!priorEvent_dist %in% unlist(ok_dists)) {
      stop("The prior distribution for the event model coefficients should be one of ",
           paste(names(ok_dists), collapse = ", "))
    } else if (priorEvent_dist %in% c("normal", "t")) {
      priorEvent_dist <- ifelse(priorEvent_dist == "normal", 1L, 2L)
      # !!! Need to think about whether 2.5 is appropriate default value here
      priorEvent_scale <- rstanarm:::set_prior_scale(priorEvent_scale, default = 2.5, 
                                     link = "none")
    } else {
      priorEvent_dist <- ifelse(priorEvent_dist == "hs", 3L, 4L)
    }
    
    priorEvent_df <- rstanarm:::maybe_broadcast(priorEvent_df, e_K)
    priorEvent_df <- as.array(pmin(.Machine$double.xmax, priorEvent_df))
    priorEvent_mean <- rstanarm:::maybe_broadcast(priorEvent_mean, e_K)
    priorEvent_mean <- as.array(priorEvent_mean)
    priorEvent_scale <- rstanarm:::maybe_broadcast(priorEvent_scale, e_K)
    priorEvent_scale <- as.array(priorEvent_scale)
  }
  if (is.null(priorEvent_intercept)) {
    priorEvent_dist_for_intercept <- 0L
    priorEvent_mean_for_intercept <- 0 
    priorEvent_scale_for_intercept <- priorEvent_df_for_intercept <- 1
  } else {
    if (!is.list(priorEvent_intercept)) 
      stop("'priorEvent_intercept' should be a named list.")
    priorEvent_dist_for_intercept <- priorEvent_intercept$dist
    priorEvent_scale_for_intercept <- priorEvent_intercept$scale
    priorEvent_mean_for_intercept <- priorEvent_intercept$location
    priorEvent_df_for_intercept <- priorEvent_intercept$df 
    priorEvent_df_for_intercept[is.na(priorEvent_df_for_intercept)] <- 1
    
    if (!priorEvent_dist_for_intercept %in% unlist(ok_intercept_dists))
      stop("The prior distribution for the event model intercept should be one of ",
           paste(names(ok_intercept_dists), collapse = ", "))
    priorEvent_dist_for_intercept <- 
      ifelse(priorEvent_dist_for_intercept == "normal", 1L, 2L)
    # !!! Need to think about whether 10 is appropriate default value here      
    priorEvent_scale_for_intercept <- 
      rstanarm:::set_prior_scale(priorEvent_scale_for_intercept, default = 10, 
                      link = "none")
    priorEvent_df_for_intercept <- min(.Machine$double.xmax, priorEvent_df_for_intercept)
  }

  # Priors for association parameters
  priorAssoc_scaled <- priorAssoc_ops$scaled
  priorAssoc_min_prior_scale <- priorAssoc_ops$min_prior_scale
  
  if (is.null(priorAssoc)) {
    priorAssoc_dist <- 0L
    priorAssoc_mean <- as.array(rep(0, a_K))
    priorAssoc_scale <- priorAssoc_df <- as.array(rep(1, a_K))
  } else {
    if (!is.list(priorAssoc)) 
      stop("'priorAssoc' should be a named list.")
    priorAssoc_dist <- priorAssoc$dist
    priorAssoc_scale <- priorAssoc$scale
    priorAssoc_mean <- priorAssoc$location
    priorAssoc_df <- priorAssoc$df
    priorAssoc_df[is.na(priorAssoc_df)] <- 1
    if (!priorAssoc_dist %in% unlist(ok_dists)) {
      stop("The prior distribution for the event model coefficients should be one of ",
           paste(names(ok_dists), collapse = ", "))
    } else if (priorAssoc_dist %in% c("normal", "t")) {
      priorAssoc_dist <- ifelse(priorAssoc_dist == "normal", 1L, 2L)
      # !!! Need to think about whether 2.5 is appropriate default value here
      priorAssoc_scale <- rstanarm:::set_prior_scale(priorAssoc_scale, default = 2.5, 
                                     link = "none")      
    } else {
      priorAssoc_dist <- ifelse(priorAssoc_dist == "hs", 3L, 4L)
    }
    
    priorAssoc_df <- rstanarm:::maybe_broadcast(priorAssoc_df, a_K)
    priorAssoc_df <- as.array(pmin(.Machine$double.xmax, priorAssoc_df))
    priorAssoc_mean <- rstanarm:::maybe_broadcast(priorAssoc_mean, a_K)
    priorAssoc_mean <- as.array(priorAssoc_mean)
    priorAssoc_scale <- rstanarm:::maybe_broadcast(priorAssoc_scale, a_K)
    priorAssoc_scale <- as.array(priorAssoc_scale)
  }
  
  # Minimum scaling of priors for longitudinal submodel(s)
  if (priorLong_scaled && priorLong_dist > 0L) {
    for (m in 1:M) {
      if (y_K[m] > 0L) {
        if (m == 1L) {
          mark_start <- 1 
          mark_end <- y_K[1]
        } else {
          mark_start <- sum(y_K[1:(m-1)]) + 1
          mark_end <- sum(y_K[1:m])
        }
        if (is_gaussian[[m]]) {
          ss <- 2 * sd(y[[m]])
          priorLong_scale[mark_start:mark_end] <- ss * priorLong_scale[mark_start:mark_end]
          priorLong_scale_for_intercept[[m]] <-  ss * priorLong_scale_for_intercept[[m]]
        }
        if (!QR) 
          priorLong_scale[mark_start:mark_end] <- 
            pmax(priorLong_min_prior_scale, priorLong_scale[mark_start:mark_end] / 
                 apply(xtemp[[m]], 2L, FUN = function(x) {
                   num.categories <- length(unique(x))
                   x.scale <- 1
                   if (num.categories == 2) x.scale <- diff(range(x))
                   else if (num.categories > 2) x.scale <- 2 * sd(x)
                   return(x.scale)
                 }))      
      }

    }
  }
  priorLong_scale <- as.array(pmin(.Machine$double.xmax, priorLong_scale))
  priorLong_scale_for_intercept <- 
    as.array(pmin(.Machine$double.xmax, priorLong_scale_for_intercept))

  # Minimum scaling of priors for event submodel
  if (priorEvent_scaled && priorEvent_dist > 0L) {
    priorEvent_scale <- pmax(priorEvent_min_prior_scale, priorEvent_scale / 
                            apply(e_xtemp, 2L, FUN = function(x) {
                              num.categories <- length(unique(x))
                              e.x.scale <- 1
                              if (num.categories == 2) e.x.scale <- diff(range(x))
                              else if (num.categories > 2) e.x.scale <- 2 * sd(x)
                              return(e.x.scale)
                            }))
  }
  priorEvent_scale <- as.array(pmin(.Machine$double.xmax, priorEvent_scale))
  priorEvent_scale_for_intercept <- 
    min(.Machine$double.xmax, priorEvent_scale_for_intercept)

  # Minimum scaling of priors for association parameters    
  if (priorAssoc_dist > 0L) {
    priorAssoc_scale <- pmax(priorAssoc_min_prior_scale, priorAssoc_scale)
  }
  priorAssoc_scale <- as.array(pmin(.Machine$double.xmax, priorAssoc_scale))

  # QR not yet implemented for stan_jm  
  if (QR) {
    stop("QR decomposition is not yet supported by stan_jm")
    if (ncol(xtemp) <= 1)
      stop("'QR' can only be specified when there are multiple predictors.")
    cn <- colnames(xtemp)
    decomposition <- qr(xtemp)
    sqrt_nm1 <- sqrt(nrow(xtemp) - 1L)
    Q <- qr.Q(decomposition)
    R_inv <- qr.solve(decomposition, Q) * sqrt_nm1
    xtemp <- Q * sqrt_nm1
    colnames(xtemp) <- cn
    xbar <- c(xbar %*% R_inv)
  }
 
  #=========================
  # Data for export to Stan
  #=========================

  standata <- list(  
    # dimensions
    M = as.integer(M),
    Npat = as.integer(Npat),
    y_N = as.array(y_N), 
    y_real_N = as.array(y_real_N), 
    y_int_N = as.array(y_int_N), 
    y_N01 = as.array(t(sapply(y_N01, cbind))), 
    y_K = as.array(y_K), 
    sum_y_N = as.integer(sum_y_N),
    sum_y_real_N = as.integer(sum_y_real_N),
    sum_y_int_N = as.integer(sum_y_int_N),
    sum_y_K = as.integer(sum_y_K),
    e_K = as.integer(e_K),
    a_K = as.integer(a_K),
    quadnodes = quadnodes,
    Npat_times_quadnodes = as.integer(Npat * quadnodes),
    sum_y_has_intercept = as.integer(sum_y_has_intercept), 
    sum_y_has_intercept_unbound = as.integer(sum_y_has_intercept_unbound), 
    sum_y_has_intercept_lobound = as.integer(sum_y_has_intercept_lobound), 
    sum_y_has_intercept_upbound = as.integer(sum_y_has_intercept_upbound), 
    sum_y_has_dispersion = as.integer(sum_y_has_dispersion),
    
    # data for longitudinal submodel(s)
    link = as.array(link),
    y_centre = as.integer(centreLong),
    y_has_intercept = as.array(as.integer(y_has_intercept)),
    y_has_intercept_unbound = as.array(y_has_intercept_unbound),
    y_has_intercept_lobound = as.array(y_has_intercept_lobound),
    y_has_intercept_upbound = as.array(y_has_intercept_upbound),
    y_has_weights = as.array(y_has_weights),
    y_has_dispersion = as.array(as.numeric(y_has_dispersion)),
    y_real = as.array(as.numeric(unlist(y[y_is_real]))),
    y_int = as.array(as.integer(unlist(y[!y_is_real]))),
    y_beg = as.array(y_beg),  # indexing for combined response vector
    y_end = as.array(y_end),
    y_real_beg = as.array(y_real_beg),
    y_real_end = as.array(y_real_end),
    y_int_beg = as.array(y_int_beg), 
    y_int_end = as.array(y_int_end),
    y_weights = as.array(do.call(c, y_weights)),
    trials = as.array(as.integer(do.call(c, trials))),
    y_xbar = if (centreLong) as.array(do.call(c, xbar)) else double(0),
    y_X = as.array(as.matrix(Matrix::bdiag(xtemp))),
    
    # data for event submodel
    basehaz_weibull = as.integer(base_haz_weibull),
    basehaz_piecewise = as.integer(base_haz_piecewise),
    basehaz_splines = as.integer(base_haz_splines),
    splines_df = as.integer(splines_df),
    e_centre = as.integer(centreEvent),
    e_has_intercept = as.integer(e_has_intercept),
    nrow_y_Xq = NROW(xqtemp[[1]]),
    nrow_e_Xq = NROW(e_xtemp),
    y_Xq = as.array(as.matrix(Matrix::bdiag(xqtemp))),
    e_Xq = e_xtemp,
    e_times = c(eventtime, unlist(quadpoint)),
    e_d = c(d, rep(1, length(unlist(quadpoint)))),
    e_xbar = if (centreEvent) as.array(e_xbar) else double(0),
    quadweight_times_half_eventtime = weights_times_half_eventtime,
    
    # data for association structure
    assoc = as.integer(a_K > 0L),
    has_assoc_ev = as.array(as.integer(has_assoc$etavalue)),
    has_assoc_es = as.array(as.integer(has_assoc$etaslope)),
    has_assoc_cv = as.array(as.integer(has_assoc$muvalue)),
    has_assoc_cs = as.array(as.integer(has_assoc$muslope)),
    sum_has_assoc_ev = as.integer(sum(has_assoc$etavalue)),
    sum_has_assoc_es = as.integer(sum(has_assoc$etaslope)),
    sum_has_assoc_cv = as.integer(sum(has_assoc$muvalue)),
    sum_has_assoc_cs = as.integer(sum(has_assoc$muslope)),
    sum_size_which_b = as.integer(sum(size_which_b)),
    size_which_b = as.array(size_which_b),
    which_b = as.array(unlist(which_b)),

    # priors
    priorLong_dist = priorLong_dist, 
    priorLong_dist_for_intercept = priorLong_dist_for_intercept,  
    priorEvent_dist = priorEvent_dist,
    priorEvent_dist_for_intercept = priorEvent_dist_for_intercept,
    priorAssoc_dist = priorAssoc_dist,    
    
    # hyperparameters for priors
    priorLong_mean = priorLong_mean, 
    priorLong_mean_for_intercept = priorLong_mean_for_intercept,
    priorEvent_mean = priorEvent_mean, 
    priorEvent_mean_for_intercept = priorEvent_mean_for_intercept,
    priorAssoc_mean = priorAssoc_mean, 
    priorLong_scale = priorLong_scale, 
    priorLong_scale_for_intercept = priorLong_scale_for_intercept, 
    priorEvent_scale = priorEvent_scale, 
    priorEvent_scale_for_intercept = priorEvent_scale_for_intercept, 
    priorAssoc_scale = priorAssoc_scale, 
    priorLong_df = priorLong_df, 
    priorLong_df_for_intercept = priorLong_df_for_intercept,  
    priorEvent_df = priorEvent_df, 
    priorEvent_df_for_intercept = priorEvent_df_for_intercept,
    priorAssoc_df = priorAssoc_df, 
    priorLong_scale_for_dispersion = as.array(priorLong_scale_for_dispersion),
    priorEvent_scale_for_weibull = 
      if (base_haz_weibull) as.array(priorEvent_scale_for_weibull) else as.array(double(0)),
    
    prior_PD = as.integer(prior_PD)
  )  
  
  # data for random effects
  group <- lapply(1:M, function(x) {
                    pad_reTrms(Z = Z[[x]], 
                               cnms = y_cnms[[x]], 
                               flist = y_flist[[x]])})
  Z     <- lapply(1:M, function(x) group[[x]]$Z)
  y_cnms <- lapply(1:M, function(x) group[[x]]$cnms)
  y_flist_padded <- lapply(1:M, function(x) group[[x]]$flist)
  t <- length(cnms_nms) # num. of unique grouping factors
  t_i <- which(cnms_nms == id_var) # index of patient-level grouping factor
  pmat <- matrix(0, t, M)
  for (i in 1:t) {
    for (j in 1:M) {
      pmat[i,j] <- length(y_cnms[[j]][[cnms_nms[i]]])
    }
  }
  p <- rowSums(pmat)
  lmat <- matrix(0, t, M)
  l <- c()
  for (i in 1:t) {
    for (j in 1:M) {
      lmat[i,j] <- nlevels(y_flist_padded[[j]][[cnms_nms[i]]])
    }
    l[i] <- max(lmat[i,])
    if (!all(lmat[i,] %in% c(0, l[i])))
      stop("The number of factor levels for each of the grouping factors ",
           "must be the same in each of the longitudinal submodels")
  }
  qmat <- l * pmat
  q <- l * p
  
  # Names of clustering variables
  group_nms <- lapply(y_cnms, names)
  # Names of random effects and random coefficients
  b_nms <- character()
  g_nms <- character() 
  for (m in 1:M) {
    for (i in seq_along(group_nms[[m]])) {
      # !!! if you change this change .pp_data_mer_z() as well
      nm <- group_nms[[m]][i]
      nms_i <- paste(y_cnms[[m]][[nm]], nm)
      if (length(nms_i) == 1) {
        b_nms <- c(b_nms, paste0("Long", m, "|", nms_i, ":", levels(y_flist_padded[[m]][[nm]])))
      } else {
        b_nms <- c(b_nms, c(t(sapply(nms_i, function(x) 
                          paste0("Long", m, "|", x, ":", levels(y_flist_padded[[m]][[nm]]))))))
      }
      g_nms <- c(g_nms, paste0("Long", m, "|", nms_i)) 
    }
  }
  standata$t <- t
  standata$t_i <- as.integer(t_i)
  standata$pmat <- as.array(pmat)
  standata$p <- as.array(p)
  standata$l <- as.array(l)
  standata$qmat <- as.array(qmat)
  standata$q <- as.array(q)
  standata$len_theta_L <- sum(choose(p, 2), p)
  Zmerge <- Matrix::bdiag(Z)
  standata$len_b <- as.integer(ncol(Zmerge))
  parts <- rstan::extract_sparse_parts(Zmerge)
  standata$num_non_zero <- as.integer(length(parts$w))
  standata$w <- parts$w
  standata$v <- parts$v
  standata$u <- parts$u
 
  # data for random effects in GK quadrature
  Zqmerge <- Matrix::bdiag(Zq)
  parts_Zq <- rstan::extract_sparse_parts(Zqmerge)
  standata$num_non_zero_Zq <- as.integer(length(parts_Zq$w))
  standata$w_Zq <- parts_Zq$w
  standata$v_Zq <- parts_Zq$v
  standata$u_Zq <- parts_Zq$u

  # data for calculating eta slope in GK quadrature 
  standata$y_dXdtq = as.array(as.matrix(Matrix::bdiag(dxdtqtemp)))
  dZdtqmerge <- Matrix::bdiag(dZdtq)
  parts_dZdtq <- rstan::extract_sparse_parts(dZdtqmerge)
  standata$num_non_zero_dZdtq <- as.integer(length(parts_dZdtq$w))
  standata$w_dZdtq <- parts_dZdtq$w
  standata$v_dZdtq <- parts_dZdtq$v
  standata$u_dZdtq <- parts_dZdtq$u    
  
  # hyperparameters for random effects model
  decov_args <- prior_covariance
  standata$shape <- as.array(rstanarm:::maybe_broadcast(decov_args$shape, t))
  standata$scale <- as.array(rstanarm:::maybe_broadcast(decov_args$scale, t))
  standata$len_concentration <- sum(p[p > 1])
  standata$concentration <- 
    as.array(rstanarm:::maybe_broadcast(decov_args$concentration, sum(p[p > 1])))
  standata$len_regularization <- sum(p > 1)
  standata$regularization <- 
    as.array(rstanarm:::maybe_broadcast(decov_args$regularization, sum(p > 1))) 

  standata$family <- as.array(sapply(1:M, function(x) {
                       return_fam <- switch(family[[x]]$family, 
                            gaussian = 1L, 
                            Gamma = 2L,
                            inverse.gaussian = 3L,
                            binomial = 4L,
                            poisson = 6L,
                            "neg_binomial_2" = 7L)
                       if (is_bernoulli[[x]]) return_fam <- 5L
                       return_fam}))
  standata$any_fam_3 <- as.integer(any(standata$family == 3L))
  
  standata$e_ns_times <- if (base_haz_splines) 
    as.array(splines::ns(standata$e_times, splines_df)) else as.array(matrix(0,0,0))
  
  #================
  # Initial values
  #================
 
  if (init == "model_based") {
    y_gamma_unbound <- unlist(y_gamma_unbound)
    y_gamma_lobound <- unlist(y_gamma_lobound)
    y_gamma_upbound <- unlist(y_gamma_upbound)
    y_z_beta <- (unlist(y_beta) - priorLong_mean) / priorLong_scale
    y_dispersion_unscaled <- y_dispersion / priorLong_scale_for_dispersion
    e_z_beta <- (e_beta - priorEvent_mean) / priorEvent_scale 

    y_hs <- if (priorLong_dist <= 2L) 0 
            else if (priorLong_dist == 3L) 2
            else if (priorLong_dist == 4L) 4
    e_hs <- if (priorEvent_dist <= 2L) 0 
            else if (priorEvent_dist == 3L) 2
            else if (priorEvent_dist == 4L) 4
    a_hs <- if (priorAssoc_dist <= 2L) 0 
            else if (priorAssoc_dist == 3L) 2
            else if (priorAssoc_dist == 4L) 4

    if (prior_covariance$dist == "decov") {
      len_z_T <- 0
      for (i in 1:t) {
        if (p[i] > 2) 
          for (j in 3:p[i]) len_z_T <- len_z_T + p[i] - 1;
      }
    } else if (prior_covariance$dist == "lkjcorr") { 
      sd_b_unscaled <- (unlist(sd_b)) / prior_scale_for_sd_b
      
    }
    
    model_based_inits <- Filter(function(x) (!is.null(x)), c(
      list(
      y_gamma_unbound = if (sum_y_has_intercept_unbound) as.array(y_gamma_unbound) else double(0),
      y_gamma_lobound = if (sum_y_has_intercept_lobound) as.array(y_gamma_lobound) else double(0),
      y_gamma_upbound = if (sum_y_has_intercept_upbound) as.array(y_gamma_upbound) else double(0),
      y_z_beta = if (sum_y_K) as.array(y_z_beta) else double(0),
      y_dispersion_unscaled = if (sum_y_has_dispersion) as.array(y_dispersion_unscaled) else double(0),
      e_gamma = if (e_has_intercept) as.array(0) else double(0),
      splines_coefs = if (!is.null(df)) as.array(rep(0, splines_df)) else double(0),
      e_z_beta = if (e_K) as.array(e_z_beta) else double(0),
      weibull_shape_unscaled = if (base_haz_weibull) 
        as.array(runif(1, 0.5, 3) / priorEvent_scale_for_weibull) else double(0),
      a_z_beta = if (a_K) as.array(rep(0, a_K)) else double(0),
      z_b = as.array(runif(standata$len_b, -0.5, 0.5)),
      y_global = as.array(runif(y_hs)),
      y_local = if (y_hs) as.array(runif(y_hs * sum_y_K), dim = c(y_hs, sum_y_K))
        else matrix(0,0,0),
      e_global = as.array(runif(e_hs)),
      e_local = if (e_hs) as.array(runif(e_hs * e_K), dim = c(e_hs, e_K))
        else matrix(0,0,0),
      a_global = as.array(runif(a_hs)),
      a_local = if (a_hs) as.array(runif(a_hs * a_K), dim = c(a_hs, a_K))
        else matrix(0,0,0)      
      ),
      if (prior_covariance$dist == "decov") list(
      z_T = as.array(runif(len_z_T, -0.5, 0.5)),
      rho = as.array(runif(sum(p)-t)),
      zeta = as.array(runif(standata$len_concentration)),
      tau = as.array(runif(t))
      ),
      if (prior_covariance$dist == "lkjcorr") list(
      sd_b_unscaled = as.array(sd_b_unscaled),
      L_b_Corr = t(chol(as.matrix(Matrix::bdiag(b_Corr))))
      )
    ))
    init <- function() model_based_inits
  }


  #===========
  # Fit model
  #===========

  # !!! REMOVE
  standata <<- standata

  # call stan() to draw from posterior distribution
  stanfit <- stanmodels$jm
  pars <- c(if (sum_y_has_intercept_unbound) "y_gamma_unbound",
            if (sum_y_has_intercept_lobound) "y_gamma_lobound",
            if (sum_y_has_intercept_upbound) "y_gamma_upbound",
            #if (sum_y_has_intercept_unbound) "y_alpha_unbound", 
            #if (sum_y_has_intercept_lobound) "y_alpha_lobound", 
            #if (sum_y_has_intercept_upbound) "y_alpha_upbound", 
            if (sum_y_K) "y_beta",
            if (e_has_intercept) "e_gamma",
            if (e_K) "e_beta",
            if (a_K) "a_beta",
            if (Npat) "b_by_model",
            "y_dispersion", 
            if (base_haz_weibull) "weibull_shape",
            if (base_haz_splines) "splines_coefs")
           #"mean_PPD"

  #cat("\n--> Fitting joint model now...")
  cat("\nPlease note the warmup phase may be much slower than",
      "later iterations!\n")             
  sampling_args <- rstanarm:::set_sampling_args(
    object = stanfit, 
    prior = priorLong, # determines default adapt_delta value?
    user_dots = list(...), 
    user_adapt_delta = adapt_delta, 
    data = standata, 
    pars = pars, 
    init = init,
    show_messages = FALSE)
  stanfit <- do.call(sampling, sampling_args)

  if (QR) {  # not yet implemented for stan_jm
    thetas <- extract(stanfit, pars = "beta", inc_warmup = TRUE, 
                      permuted = FALSE)
    betas <- apply(thetas, 1:2, FUN = function(theta) R_inv %*% theta)
    end <- tail(dim(betas), 1L)
    for (chain in 1:end) for (param in 1:nrow(betas)) {
      stanfit@sim$samples[[chain]][[has_intercept + param]] <-
        if (ncol(xtemp) > 1) betas[param, , chain] else betas[param, chain]
    }
  }

  # Names for coefs from submodel(s)
  int_nms <- unlist(lapply(1:M, function(x) 
                    if (y_has_intercept[x]) paste0("Long", x, "|(Intercept)")))
  y_nms   <- unlist(lapply(1:M, function(x) 
                    paste0("Long", x, "|", colnames(xtemp[[x]]))))
  e_nms   <- paste0("Event|", colnames(e_x))    
  
  # Names for vector of association parameters
  a_nms <- character()  
  for (m in 1:M) {
    if (has_assoc$etavalue[m]) a_nms <- c(a_nms, paste0("Assoc|Long", m,":eta value"))
    if (has_assoc$etaslope[m]) a_nms <- c(a_nms, paste0("Assoc|Long", m,":eta slope"))
    if (has_assoc$muvalue[m]) a_nms <- c(a_nms, paste0("Assoc|Long", m,":mu value"))
    if (has_assoc$muslope[m]) a_nms <- c(a_nms, paste0("Assoc|Long", m,":mu slope"))
  }
  if (sum(size_which_b)) {
    temp_g_nms <- lapply(1:M, FUN = function(m) {
                    all_nms <- paste0(paste0("Long", m, ":b["), y_cnms[[m]][[id_var]], "]")
                    all_nms[which_b[[m]]]})
    a_nms <- c(a_nms, paste0("Assoc|", unlist(temp_g_nms)))
  }
  
  # Names for vector of dispersion parameters
  d_nms <- character()  
  for (m in 1:M) {
    if (rstanarm:::is.gaussian(famname[[m]]))   d_nms <- c(d_nms, paste0("Long", m,"|sigma"))
    else if (rstanarm:::is.gamma(famname[[m]])) d_nms <- c(d_nms, paste0("Long", m,"|shape"))
    else if (rstanarm:::is.ig(famname[[m]]))    d_nms <- c(d_nms, paste0("Long", m,"|lambda"))
    else if (rstanarm:::is.nb(famname[[m]]))    d_nms <- c(d_nms, paste0("Long", m,"|overdispersion"))
  }
                  
  new_names <- c(int_nms,
                 y_nms,
                 e_nms,
                 a_nms,                   
                 if (length(group)) c(paste0("b[", b_nms, "]")),
                 d_nms,
                 if (base_haz_weibull) "Event|weibull shape",               
                 if (base_haz_splines) paste0("Event|basehaz coef", seq(splines_df)),               
                #"mean_PPD", 
                 "log-posterior")
  stanfit@sim$fnames_oi <- new_names
  
  n_grps <- l - 1
  names(n_grps) <- cnms_nms  # n_grps is num. of levels within each grouping factor
  names(p) <- cnms_nms       # p is num. of variables within each grouping factor
  
  #colnames(Z) <- b_names(names(stanfit), value = TRUE)
  fit <- rstanarm:::nlist(stanfit, family, formula = c(formulaLong, formulaEvent), 
                          id_var, time_var, offset = NULL, base_haz, quadnodes,
                          df = if (base_haz_splines) splines_df else NULL,
                          M, cnms, y_N, y_cnms, y_flist, Npat, n_grps, assoc = has_assoc, 
                          fr = c(y_fr, list(e_fr)),
                          x = lapply(1:M, function(i) 
                            if (getRversion() < "3.2.0") cBind(x[[i]], Z[[i]]) else cbind2(x[[i]], Z[[i]])),
                          xq = lapply(1:M, function(i) 
                            if (getRversion() < "3.2.0") cBind(xq[[i]], Zq[[i]]) else cbind2(xq[[i]], Zq[[i]])),                 
                          dxdtq = lapply(1:M, function(i) 
                            if (getRversion() < "3.2.0") 
                              cBind(dxdtqtemp[[i]], dZdtq[[i]]) else cbind2(dxdtqtemp[[i]], dZdtq[[i]])),                          
                          y = y, e_x, eventtime, d,
                          standata, dataLong, dataEvent, call, terms = NULL, model = NULL,                          
                          prior.info = rstanarm:::get_prior_info(call, formals()),
                          na.action, algorithm = "sampling", init, glmod = y_mod, coxmod = e_mod)
  out <- stanjm(fit)
  
  return(out)
}


# Check the id_var argument is valid and is included appropriately in the
# formulas for each of the longitudinal submodels
#
# @param id_var The character string that the user specified for the id_var
#   argument -- will have been set to NULL if the argument was missing.
# @param y_cnms A list of length M with the cnms for each longitudinal submodel
# @return Returns the character string corresponding to the appropriate id_var.
#   This will either be the user specified id_var argument or the only grouping
#   factor.
check_id_var <- function(id_var, y_cnms) {
  len_cnms <- sapply(y_cnms, length)
  if (any(len_cnms > 1L)) {  # more than one grouping factor
    if (is.null(id_var)) {
      stop("'id_var' must be specified when using more than one grouping factor",
           call. = FALSE)
    } else {
      lapply(y_cnms, function(x)  if (!(id_var %in% names(x)))
        stop("`id_var' must be included as a grouping factor in each ",
             "of the longitudinal submodels", call. = FALSE)) 
      return(id_var)
    }      
  } else {  # only one grouping factor (assumed to be subject ID)
    only_cnm <- unique(sapply(y_cnms, names))
    if (length(only_cnm) > 1L)
      stop("The grouping factor (ie, subject ID variable) is not the ",
           "same in all longitudinal submodels", call. = FALSE)
    if ((!is.null(id_var)) && (!identical(id_var, only_cnm)))
        warning("The user specified 'id_var' (", paste(id_var), 
                ") and the assumed ID variable based on the single ",
                "grouping factor (", paste(only_cnm), ") are not the same; ", 
                "`id_var' will be ignored", call. = FALSE, immediate. = TRUE)
    return(only_cnm)
  }
}


# Check the factor list corresponding to subject ID is the same in each 
# of the longitudinal submodels
#
check_id_list <- function(id_var, y_flist) {
  id_list <- unique(lapply(y_flist, function(x) unique(x[[id_var]])))
  if (length(id_list) > 1L)
    stop("The subject IDs are not the same in all longitudinal submodels.",
         call. = FALSE)
  unlist(id_list)  
}


# Function to check that the assoc argument only includes supported association
# types. The function returns a list with logicals specifying which association
# type have been requested.
# 
# @param x The assoc argument specified by the user -- should be a character 
#   vector or NULL
# @param supported_assocs A character vector showing the supported
#   association types
# @return A list of logicals indicating the desired association types
validate_assoc <- function(x, supported_assocs) {
  
  # Identify which association types were specified
  x_tmp <- gsub("^shared_b.*", "shared_b", x)
  assoc <- sapply(supported_assocs, function(y) y %in% x_tmp, simplify = FALSE)
  if (is.null(x_tmp)) {
    assoc$null <- TRUE
  } else if (is.character(x_tmp)) {
    if (!all(x_tmp %in% supported_assocs))
      stop("An unsupported association type has been specified. The ",
           "'assoc' argument can only include the following association ", 
           "types: ", paste(supported_assocs, collapse = ", "), call. = FALSE)
    if ((assoc$null) && (length(assoc) > 1L))
      stop("In 'assoc' argument, 'null' cannot be specified in ",
           "conjuction with another association type", call. = FALSE)
    if (assoc$etavalue && assoc$muvalue)
      stop("In 'assoc' argument, 'etavalue' and 'muvalue' cannot be specified ",
           "together", call. = FALSE)
    if (assoc$etaslope && assoc$muslope)
      stop("In 'assoc' argument, 'etaslope' and 'muslope' cannot be specified ",
           "together", call. = FALSE)    
  } else { 
    stop("'assoc' argument should be a character vector or, for a multivariate ",
         "joint model, possibly a list of character vectors.", call. = FALSE)
  }
  
  # Identify which subset of shared random effects were specified
  val <- grep("^shared_b.*", x, value = TRUE)
  if (length(val)) {
    val <- unlist(strsplit(val, "shared_b"))[-1]
  }
  if (length(val)) {
    assoc$which_b <- tryCatch(eval(parse(text = paste0("c", val))), 
                            error = function(x) 
                              stop("Incorrect specification of the 'shared_b' ",
                                   "association structure. See Examples in help ",
                                   "file.", call. = FALSE))
  } else assoc$which_b <- numeric(0)

  assoc
}


# Function to return the variables used in fitting a model, as well as the
# predvars, for both the fixed and random parts of a (g)lmer fit
#
# @param mod A (g)lmer model object
# @return A named list of lists
get_formvars <- function(mod) {
  vars_f <- grep("", attr(terms(mod, fixed.only = TRUE), "variables"), value = TRUE)
  vars_r <- grep("", attr(terms(mod, random.only = TRUE), "variables"), value = TRUE)
  predvars_f <- grep("", attr(terms(mod, fixed.only = TRUE), "predvars"), value = TRUE)
  predvars_r <- grep("", attr(terms(mod, random.only = TRUE), "predvars"), value = TRUE)
  list(formvars = list(fixed = vars_f, random = vars_r),
       predvars = list(fixed = predvars_f, random = predvars_r))
}


# Function to replace the variables in the model formula with the prediction variables
#
# @param mod A (g)lmer model object from which to extract the model formula
# @param formvars A list with components fixed and random which give the original 
#   variables used in the model formula
# @param predvars A list with components fixed and random which gives the new 
#   variables used to replace the variables in formvars
# @return The reformulated model formula with the variables replaced by predvars
substitute_vars <- function(mod, formvars, predvars) {
  fm <- formula(mod)
  if (!identical(formvars, predvars)) {
    for (j in 2:length(formvars$fixed))
      fm <- gsub(formvars$fixed[[j]], predvars$fixed[[j]], fm, fixed = TRUE)    
    for (j in 2:length(formvars$random))
      fm <- gsub(formvars$random[[j]], predvars$random[[j]], fm, fixed = TRUE)    
    fm <- reformulate(fm[[3]], response = fm[[2]])
  }
  fm
}


# Function to return control arguments for lmer or glmer call
#
get_control_args <- function(glmer = FALSE, norank = FALSE) {
  if (glmer) {
    if (norank) {
      lme4::glmerControl(check.nlev.gtreq.5 = "ignore",
                         check.nlev.gtr.1 = "stop",
                         check.nobs.vs.rankZ = "ignore",
                         check.nobs.vs.nlev = "ignore",
                         check.nobs.vs.nRE = "ignore",
                         check.rankX = "ignore")      
    } else {
      lme4::glmerControl(check.nlev.gtreq.5 = "ignore",
                         check.nlev.gtr.1 = "stop",
                         check.nobs.vs.rankZ = "ignore",
                         check.nobs.vs.nlev = "ignore",
                         check.nobs.vs.nRE = "ignore")      
    }    
  } else {
    if (norank) {
      lme4::lmerControl(check.nlev.gtreq.5 = "ignore",
                        check.nlev.gtr.1 = "stop",
                        check.nobs.vs.rankZ = "ignore",
                        check.nobs.vs.nlev = "ignore",
                        check.nobs.vs.nRE = "ignore",
                        check.rankX = "ignore")       
    } else {
      lme4::lmerControl(check.nlev.gtreq.5 = "ignore",
                        check.nlev.gtr.1 = "stop",
                        check.nobs.vs.rankZ = "ignore",
                        check.nobs.vs.nlev = "ignore",
                        check.nobs.vs.nRE = "ignore")      
    }
  }
}

# Function to return standardised GK quadrature points and weights
#
# @return A list with two named elements (points and weights) each
#   of which is a numeric vector with length equal to the number of
#   quadrature nodes
get_quadpoints <- function(nodes = 15) {
  if (nodes == 15) {
    list(
      points = c(
        -0.991455371120812639207,
        -0.949107912342758524526,
        -0.86486442335976907279,
        -0.7415311855993944398639,
        -0.5860872354676911302941,
        -0.4058451513773971669066,
        -0.2077849550078984676007,
        0,
        0.2077849550078984676007,
        0.405845151377397166907,
        0.5860872354676911302941,
        0.741531185599394439864,
        0.86486442335976907279,
        0.9491079123427585245262,
        0.991455371120812639207),
      weights = c(
        0.0229353220105292249637,
        0.063092092629978553291,
        0.10479001032225018384,
        0.140653259715525918745,
        0.1690047266392679028266,
        0.1903505780647854099133,
        0.204432940075298892414,
        0.209482141084727828013,
        0.204432940075298892414,
        0.1903505780647854099133,
        0.169004726639267902827,
        0.140653259715525918745,
        0.1047900103222501838399,
        0.063092092629978553291,
        0.0229353220105292249637))      
  } else if (nodes == 11) {
    list(
      points = c(
        -0.984085360094842464496,
        -0.906179845938663992798,
        -0.754166726570849220441,
        -0.5384693101056830910363,
        -0.2796304131617831934135,
        0,
        0.2796304131617831934135,
        0.5384693101056830910363,
        0.754166726570849220441,
        0.906179845938663992798,
        0.984085360094842464496),
      weights = c(
        0.042582036751081832865,
        0.1152333166224733940246,
        0.186800796556492657468,
        0.2410403392286475866999,
        0.272849801912558922341,
        0.2829874178574912132043,
        0.272849801912558922341,
        0.241040339228647586701,
        0.186800796556492657467,
        0.115233316622473394025,
        0.042582036751081832865))     
  } else if (nodes == 7) {
    list(
      points = c(
        -0.9604912687080202834235,
        -0.7745966692414833770359,
        -0.4342437493468025580021,
        0,
        0.4342437493468025580021,
        0.7745966692414833770359,
        0.9604912687080202834235),
      weights = c(
        0.1046562260264672651938,
        0.268488089868333440729,
        0.401397414775962222905,
        0.450916538658474142345,
        0.401397414775962222905,
        0.268488089868333440729,
        0.104656226026467265194))      
  } else stop("The specified number of Gauss-Kronrod quadrature points 
              ('quadnodes') must be either 7, 11 or 15.")  
}


# Function to return the appropriate intercept type based on family and link
#
# @param family A GLM family
# @param link A link function
# @return A list specifying whether an unbounded, lower bounded, or upper 
#   bounded intercept is required
check_intercept <- function(family, link) {
  if (family == "binomial") {
    if (link == "log") {
      unbound <- 0L; lobound <- 0L; upbound <- 1L;  # binomial, log
    } else {
      unbound <- 1L; lobound <- 0L; upbound <- 0L;  # binomial, !log
    }
  } else {
    if (link == "log") {
      unbound <- 1L; lobound <- 0L; upbound <- 0L;  # gamma/inv-gaus/poisson/nb, log
    } else {
      if (family == "gaussian") {
        unbound <- 1L; lobound <- 0L; upbound <- 0L;  # gaussian, any link
      } else {
        unbound <- 0L; lobound <- 1L; upbound <- 0L;  # gamma/inv-gaus/poisson/nb, !log      
      }
    }
  }  
  rstanarm:::nlist(unbound, lobound, upbound) 
}


# Function to check if the response vector is real or integer
#
# @param family A GLM family
# @return A logical specify whether the response is real (TRUE) or integer (FALSE)
check_response_real <- function(family) {
  if (family == "binomial" || 
      family == "poisson" ||
      family == "neg_binomial_2") {
    FALSE
  } else TRUE
}


# Function to check if the submodel should include a dispersion term
#
# @param family A GLM family
# @return A logical specify whether the submodel includes a dispersion term
check_for_dispersion <- function(family) {
  if (family == "binomial" || family == "poisson") FALSE else TRUE
}


# Function to calculate the number of association parameters in the model
#
# @param has_assoc A named list specifying whether each longitudinal submodel 
#   is linked to the event outcome using each potential type of association structure
# @param which_b A list of numeric vectors indicating the random effects from each
#   longitudinal submodel that are to be used in the shared_b association structure
# @return Integer indicating the number of association parameters in the model 
get_num_assoc_pars <- function(has_assoc, which_b) {
  sel <- c("etavalue", "etaslope", "muvalue", "muslope")
  a_K <- sum(unlist(has_assoc[sel]))
  a_K <- a_K + length(unlist(which_b))
  return(a_K)
}

# Add extra level _NEW_ to each group
# 
# @param Z ranef indicator matrix
# @param cnms group$cnms
# @param flist group$flist
pad_reTrms <- function(Z, cnms, flist) {
  l <- sapply(attr(flist, "assign"), function(i) nlevels(flist[[i]]))
  p <- sapply(cnms, FUN = length)
  last <- cumsum(l * p)
  for (i in attr(flist, "assign")) {
    levels(flist[[i]]) <- c(gsub(" ", "_", levels(flist[[i]])), 
                            paste0("_NEW_", names(flist)[i]))
  }
  n <- nrow(Z)
  mark <- length(p) - 1L
  if (getRversion() < "3.2.0") {
    Z <- cBind(Z, Matrix(0, nrow = n, ncol = p[length(p)], sparse = TRUE))
    for (i in rev(head(last, -1))) {
      Z <- cBind(cBind(Z[, 1:i, drop = FALSE],
                       Matrix(0, n, p[mark], sparse = TRUE)),
                 Z[, (i+1):ncol(Z), drop = FALSE])
      mark <- mark - 1L
    }
  }
  else {
    Z <- cbind2(Z, Matrix(0, nrow = n, ncol = p[length(p)], sparse = TRUE))
    for (i in rev(head(last, -1))) {
      Z <- cbind(Z[, 1:i, drop = FALSE],
                 Matrix(0, n, p[mark], sparse = TRUE),
                 Z[, (i+1):ncol(Z), drop = FALSE])
      mark <- mark - 1L
    }
  }
  rstanarm:::nlist(Z, cnms, flist)
}

# Drop the extra reTrms from a matrix x
#
# @param x A matrix (e.g. the posterior sample or matrix of summary stats)
# @param columns Do the columns (TRUE) or rows (FALSE) correspond to the
#   variables?
unpad_reTrms <- function(x, ...) UseMethod("unpad_reTrms")
unpad_reTrms.default <- function(x, ...) {
  if (is.matrix(x))
    return(unpad_reTrms.matrix(x, ...))
  keep <- !grepl("_NEW_", names(x), fixed = TRUE)
  x[keep]
}
unpad_reTrms.matrix <- function(x, columns = TRUE, ...) {
  nms <- if (columns) 
    colnames(x) else rownames(x)
  keep <- !grepl("_NEW_", nms, fixed = TRUE)
  if (columns) x[, keep, drop = FALSE] else x[keep, , drop = FALSE]
}






#------ Functions below here are not actually used

# Function to turn the matched call for all longitudinal submodels together
# into a matched call for just one of the longitudinal submodels
#
# @param m Integer specifying which longitudinal marker/submodel to create
#   the matched call for
# @param y_mc The user inputted matched call corresponding to all longitudinal
#   submodels
# @return The matched call for fitting the longitudinal submodel for marker m
create_m_mc <- function(m, y_mc, was_list) {    
  m_mc <- y_mc
  m_mc$formula <- if (was_list$formula) y_mc$formula[[(1+m)]] else y_mc$formula
  m_mc$data <- if (was_list$data) y_mc$data[[(1+m)]] else y_mc$data
  if (!is.null(was_list$subset))   
    m_mc$subset <- if (was_list$subset) y_mc$subset[[(1+m)]] else y_mc$subset
  if (!is.null(was_list$family))   
    m_mc$family <- if (was_list$family) y_mc$family[[(1+m)]] else y_mc$family
  return(m_mc)
}
 
# Function to create cnms object which specifies the unique cnms across
# all longitudinal submodels, rather than each submodel separately
#
# @param y_cnms A list of length M with each element corresponding to the 
#   cnms object for just one of the longitudinal submodels
get_common_cnms <- function(y_cnms) {
  y_cnms_nms <- lapply(y_cnms, names)
  cnms_nms <- unique(unlist(y_cnms_nms))
  cnms <- lapply(seq_along(cnms_nms), function(i) {
    nm <- cnms_nms[i]
    unlist(lapply(1:length(y_cnms), function(m) 
      if (nm %in% y_cnms_nms[[m]]) paste0("Long", m, "|", y_cnms[[m]][[nm]])))
  })
  names(cnms) <- cnms_nms
  return(cnms)
}

# Return the number of longitudinal markers
#
# @param x A stanjm fitted model object, a formula object 
#   or list of formula objects
get_M <- function(x) {
  if (is.call(x)) x <- eval(x)
  if (is(x, "stanjm")) {
    x$n_markers
  } else if (is(x, "formula")) {
    1L
  } else if (is.list(x)) {
    if (!all(sapply(x, function(y) (is(y, "formula")))))
      stop("Argument should be a fitted model of stanjm class, ", 
           "a formula object or a list of formula objects", call. = FALSE)
    length(x)
  } else {
    stop("Argument should be a fitted model of stanjm class, ", 
         "a formula object or a list of formula objects", call. = FALSE)
  }
}

# Check the user input to the formulaLong argument
#
# @param x The formulaLong argument of the matched call
# @return A logical indicating whether the user inputted a list 
#   of formula objects
check_arg_yform <- function(x) {
  if (is.call(x)) x <- eval(x)
  if (is(x, "formula")) {
    structure(FALSE, M = 1L)
  } else if (is.list(x)) {
    if (!all(sapply(x, function(y) (is(y, "formula")))))
      stop("'formulaLong' should be a formula object or a list of formula ",
           "objects", call. = FALSE)
    structure(TRUE, M = length(x))
  } else {
    stop("'formulaLong' should be a formula object or, for a multivariate ",
         "joint model, a list of formula objects with length equal to the ",
         "desired number of longitudinal markers", call. = FALSE)
  }
}

# Check the user input to the dataLong argument
#
# @param x The dataLong argument of the matched call
# @param M The number of longitudinal markers in the joint model
# @return A logical indicating whether the user inputted a list 
#   of data frames
check_arg_ydata <- function(x, M) {
  if (is.call(x)) x <- eval(x)
  if (is.data.frame(x)) {
    FALSE
  } else if (is.list(x)) {
    if (length(x) != M)
      stop("dataLong appears to be a list of the incorrect length",
           call. = FALSE)
    TRUE
  } else {
    stop("'dataLong' should be a data frame or possibly a list of data ",
         "frames. The latter is only required when fitting a multivariate ",
         "joint model using different data for each longitudinal submodel.",
         call. = FALSE)
  }
}

# Check the user input to the subsetLong argument
#
# @param x The subsetLong argument of the matched call
# @param M The number of longitudinal markers in the joint model
# @return An integer indicating whether the user inputted a list 
#   of subsets
check_arg_ysubset <- function(x, M) {
  if (is.call(x)) x <- eval(x)
  if (is.null(x)) {
    NULL   # no subset entered
  } else if (is.list(x)) {
    if (length(x) != M)
      stop("'subsetLong' appears to be a list of the incorrect length",
           call. = FALSE)
    TRUE   # is a list
  } else if (is.vector(x)) {
    FALSE  # not a list, but valid
  } else {
    stop("'subsetLong' should be a vector or possibly a list of vectors. ",
         "The latter is only required if fitting a multivariate joint ",
         "model and using a different subset of data for each ",
         "longitudinal submodel.", call. = FALSE)
  }
}  

# Check the user input to the family argument
#
# @param x The family argument of the matched call
# @param M The number of longitudinal markers in the joint model
# @return An integer indicating whether the user inputted a list 
#   of families
check_arg_family <- function(x, M) {
  if (is.call(x)) x <- eval(x)
  if (is.null(x)) {
    NULL
  } else if (is(x, "family")) {
    FALSE  # not a list, but valid, so convert to list
  } else if (is.list(x)) {
    if (length(x) != M)
      stop("'family' should be a family function or, for a multivariate ",
           "joint model, possibly a list of family functions", call. = FALSE)
    TRUE  # is a list
  } else {
    stop("'family' should be a family function or, possibly a list of family ",
         "functions. The latter is only required when fitting a multivariate ",
         "joint model with a different family and/or link function for some ",
         "of the longitudinal submodels.", call. = FALSE)
  }
}  


# Check if object is a list, and if not, then convert it to a
# list and possibly repeat it 'len' times
#
# @param x The object to be checked and possible converted to a list
# @param len The length of the list that is to be returned
to_list <- function(x, len = 1) {
  if (is.list(x)) {
    if (length(x) == len) {
      x
    } else if (length(x) == 1L) {
      rep(x, len)
    } else {
      stop("To be broadcast, 'x' must be a list of length 1")
    }
  } else {
    rep(list(x), len)
  }
}

# Check the user input to the assoc argument and, if necessary, convert it 
# to a list of length equal to the number of longitudinal submodels
#
# @param x The user input to the assoc argument
# @param M Integer specifying the number of longitudinal markers
# @return A list with the assoc argument for each longitudinal submodel
check_arg_assoc <- function(x, M) { 
  if (is.list(x)) {  # if list, then check length
    if (!(length(x) %in% c(1,M)))
      stop("`assoc' should be a list of length 1 or length equal to the ",
           "number of longitudinal markers", call. = FALSE)
  } else if (is.character(x)) {  # if not list, then convert to list
    x <- list(x)
  } else {  # else return error
    stop("'assoc' should be a character string or character vector or, for a ",
         "multivariate joint model, possibly a list of character strings ",
         "or character vectors. The latter is only required if using a different ",
         "association structure for linking each longitudinal submodel to the ",
         "event outcome.", call. = FALSE)
  }
  if (length(x) != M) x <- rep(x, M)
  x
} 
