# Part of the rstanjm package
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
#
#' Simulated datasets used for \pkg{rstanjm} tests
#' 
#' Simulated datasets used for running tests on the \pkg{rstanjm} package.
#'
#' @name simulated-datasets
#' @keywords internal
#' @aliases simdata_jm_cvassoc simdata_jm_cvassoc_id simdata_mjm_cvassoc
#'   simdata_mjm_cvassoc_id
#'      
#' @format 
#' \describe{
#' The simulated datasets each include 250 individuals. A binary
#' treatment code is assigned to each individual, with 
#' probability of 0.5. 
#' 
#' The datasets with \code{jm} in the name use a 
#' univariate joint model in the data generating process,
#' whilst the datasets with \code{mjm} in the name use
#' a bivariate joint model for the data generating process.   
#'   
#' The datasets without the \code{_id} suffix contain both  
#' longitudinal and event time data (they include one or more 
#' rows per individual). The datasets with the \code{_id} suffix 
#' contain event time data only (they include only one row per 
#' individual). \cr
#' \cr
#' The variable definitions are as follows: \cr
#'   \item{\code{id}}{numeric ID unique to each individual}
#'   \item{\code{trt}}{binary treatment variable taking the 
#'     values \code{0} or \code{1}}   
#'   \item{\code{t0}}{start of the time interval, taken to 
#'     be the measurement time for the longitudinal marker(s)
#'     \code{y1}, \code{y2}, etc.}
#'   \item{\code{t}}{end of the time interval, taken to be the 
#'     event time for the event indicator \code{d}}
#'   \item{\code{y1}, \code{y2}, ...}{the longitudinal 
#'     response of type gaussian (\code{_gaus}), poisson 
#'     (\code{_pois}), or bernoulli (\code{_bern}).}
#'   \item{\code{d}}{the event indicator}
#' }
#' 
NULL
