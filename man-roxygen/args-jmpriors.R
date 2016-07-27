#' @param priorLong,priorEvent,priorAssoc The prior distributions for the 
#'   regression coefficients in the longitudinal submodel(s), event submodel,
#'   and the association parameter(s). 
#'   Can be a call to \code{normal}, \code{student_t},
#'   \code{cauchy}, \code{hs} or \code{hs_plus}. See \code{\link{priors}} for
#'   details. To to omit a prior ---i.e., to use a flat (improper) uniform
#'   prior--- set equal to \code{NULL}.
#' @param priorLong_intercept,priorEvent_intercept The prior distributions  
#'   for the intercepts in the longitudinal submodel(s) and event submodel.
#'   Can be a call to \code{normal}, \code{student_t} or
#'   \code{cauchy}. See \code{\link{priors}} for details. To to omit a prior
#'   ---i.e., to use a flat (improper) uniform prior--- set
#'   equal to \code{NULL}. (\strong{Note:} If \code{centreLong} or 
#'   \code{centreEvent} is set to \code{TRUE} then the prior
#'   distribution for the intercept is set so it applies to the value when all
#'   predictors are centered.)
#' @param priorLong_ops Additional options related to prior distributions for
#'   the longitudinal submodel(s). Can  be \code{NULL} to omit a prior on the 
#'   dispersion parameters in the longitudinal submodel(s) and see 
#'   \code{\link{priorLong_options}} otherwise.
#' @param priorEvent_ops Additional options related to prior distributions 
#'   for the event submodel. Can be \code{NULL} to omit a prior on the 
#'   Weibull shape parameter (if a Weibull baseline hazard is used) and 
#'   see \code{\link{priorEvent_options}} otherwise. 
#' @param priorAssoc_ops Additional options related to prior distributions 
#'   for the association parameter(s). See \code{\link{priorAssoc_options}}.
