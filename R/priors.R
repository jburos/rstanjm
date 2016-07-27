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

#' Additional options related to prior distributions for fitting Bayesian 
#' joint models via Stan
#' 
#' These functions are used to specify additional options related to
#' prior distributions when fitting a joint model using the \code{link{stan_jm}} 
#' function, which is the main modeling function in the \pkg{rstanjm} package. 
#' However, the main functions used for specifying prior distributions for 
#' parameters in the joint model are called from the \pkg{rstanarm} package; 
#' see \code{link[rstanarm]{priors}} for further details.
#' 
#' @export 
#' @name priors
#'   
#' @details The functions shown here are used to generate additional options
#'   related to prior distributions for parameters in a joint longitudinal 
#'   and time-to-event model, estimated using the \code{link{stan_jm}}
#'   modelling function. The functions are used to specify options separately
#'   for the longitudinal submodel(s), the event submodel, and the association
#'   parameters.
#'   
#'   The prior distributions for the regression coefficients in each of the 
#'   submodels as well as the prior distribution on the covariance matrix for
#'   the group-specific parameters (or subject-specific random effects) are
#'   implemented via the \pkg{rstanarm} package; see \code{link[rstanarm]{priors}} 
#'   for details.    
#' @return A named list to be used internally by the \code{link{stan_jm}}
#'   modelling function.
#' 

#' @rdname priors
#' @export 
#' @param prior_scale_for_dispersion Prior scale for the standard error of the 
#'   regression in Gaussian models, which is given a half-Cauchy prior truncated
#'   at zero.
#' @param min_prior_scale Minimum prior scale for the intercept and 
#'   coefficients.
#' @param scaled A logical scalar, defaulting to \code{TRUE}. If \code{TRUE} the
#'   \code{prior_scale} is further scaled by the range of the predictor if the 
#'   predictor has exactly two unique values and scaled by twice the standard
#'   deviation of the predictor if it has more than two unique values.
#'
priorLong_options <- function(prior_scale_for_dispersion = 5, 
                              min_prior_scale = 1e-12, 
                              scaled = TRUE) {
  rstanarm:::validate_parameter_value(prior_scale_for_dispersion)
  rstanarm:::validate_parameter_value(min_prior_scale)
  rstanarm:::nlist(scaled, min_prior_scale, prior_scale_for_dispersion)
}

#' @rdname priors
#' @export 
#' @param prior_scale_for_weibull Prior scale for the shape parameter of the 
#'   Weibull distribution for the baseline hazard in \code{\link{stan_glm}}, 
#'   which is given a half-Cauchy truncated at zero.
#'
priorEvent_options <- function(prior_scale_for_weibull = 5,
                               min_prior_scale = 1e-12, 
                               scaled = TRUE) {
  rstanarm:::validate_parameter_value(prior_scale_for_weibull)
  rstanarm:::validate_parameter_value(min_prior_scale)
  rstanarm:::nlist(scaled, min_prior_scale, prior_scale_for_weibull)
}

#' @rdname priors
#' @export 
#'
priorAssoc_options <- function(min_prior_scale = 1e-12) {
  rstanarm:::validate_parameter_value(min_prior_scale)
  rstanarm:::nlist(min_prior_scale)
}

