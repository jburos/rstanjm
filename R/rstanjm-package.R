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

#' Bayesian joint longitudinal and time-to-event models via Stan
#'
#' @docType package
#' @name rstanjm-package
#' @aliases rstanjm
#' @useDynLib rstanjm, .registration = TRUE 
#' 
#' @import methods
#' @importFrom rstan optimizing sampling vb constrain_pars extract
#'   extract_sparse_parts get_posterior_mean stanc
#' @import stats
#' @import Rcpp
#' @description This package allows the user to fit joint models for longitudinal 
#' and time-to-event data (also known as shared parameter models). Both 
#' univariate (one longitudinal marker) and multivariate (more than one longitudinal 
#' marker) joint models are allowed. Both continuous and non-continuous (e.g.
#' binary or count data) longitudinal markers can be accomodated through a range
#' of possible link functions and error distributions. Multi-level clustered data 
#' can be accomodated (for example, patients within clinics). Competing risk events 
#' are not currenly allowed. \cr
#' \cr  
#' For each longitudinal submodel, a generalised linear mixed model is assumed with 
#' one of the family and link combinations allowed by the \pkg{\link{lme4}} package.
#' If a multivariate joint model is specified (i.e. more than one longitudinal marker)
#' then dependence between the different longitudinal markers is captured through 
#' correlated random effects; that is, a multivariate normal distribution with an 
#' unstructured correlation matrix is assumed for all random effects within a given 
#' clustering level (e.g. individual-level random effects). \cr
#' \cr
#' For the time-to-event submodel the only option currently available is a Weibull 
#' proportional hazards model. Time-varying covariates are allowed in both the 
#' longitudinal and event submodels. Competing risk events are not 
#' currently allowed. The association structure for the joint model 
#' can be based on any of the following parameterisations: 
#'    current expected value in the longitudinal submodel(s);
#'    current value of the linear predictor in the longitudinal submodel(s);
#'    first derivative (slope) for the linear predictor in the longitudinal submodel(s); 
#'    shared random effects;
#'    no association structure (equivalent to fitting separate longitudinal and survival models). \cr
#' \cr
#' Estimation of the joint model is under a Bayesian framework using Markov chain 
#' Monte Carlo (MCMC) methods, with the underlying estimation implemented using 
#' Stan \url{http://mc-stan.org/}. A variety of prior distributions are available
#' for the regression coefficients, including normal, Cauchy and t distributions, 
#' as well as shrinkage priors. The prior distribution for the covariance matrix 
#' of the individual-level (or other clustering-level) random effects is based on 
#' a novel decomposition. All prior distributions are implemented via the 
#' \pkg{rstanarm} package; see \code{\link[rstanarm]{priors}} for details.
#' 
#' @seealso \code{\link{stanjm-object}} and \code{\link{stanjm-methods}} 
#'   for details on the fitted model object returned by the \code{\link{stan_jm}}
#'   modelling function..
#'   
#'   \url{http://mc-stan.org/} for more information on the Stan C++ package used
#'   by \pkg{rstanjm} for model fitting.
#'   
#' @template reference-stan-manual
#'   
NULL
