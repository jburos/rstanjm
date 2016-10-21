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

#' Posterior uncertainty intervals
#' 
#' The \code{posterior_interval} function computes Bayesian posterior  
#' uncertainty intervals for the model parameters, via a call to the 
#' \code{posterior_interval} function in the \pkg{rstanarm} package.
#' See \code{\link[rstanarm]{posterior_interval}} for details.
#' 
#' @export
#' @templateVar stanjmArg object
#' @template args-stanjm-object
#' @param ... Arguments passed to the \code{\link[rstanarm]{posterior_interval}} 
#'   function of the \pkg{rstanarm} package.
#' 
#' @return A matrix with two columns and as many rows as model parameters (or 
#'   the subset of parameters specified by \code{pars} and/or 
#'   \code{regex_pars}). For a given value of \code{prob}, \eqn{p}, the columns 
#'   correspond to the lower and upper \eqn{100p}\% interval limits and have the
#'   names \eqn{100\alpha/2}\% and \eqn{100(1 - \alpha/2)}\%, where \eqn{\alpha 
#'   = 1-p}. For example, if \code{prob=0.9} is specified (a \eqn{90}\%
#'   interval), then the column names will be \code{"5\%"} and \code{"95\%"},
#'   respectively.
#' 
#' @template reference-gelman-carlin
#' @template reference-morey
#' 
#' @examples 
#' if (!exists("example_jm")) example(example_jm)
#' posterior_interval(example_jm)
#' posterior_interval(example_jm, regex_pars = "Long1")
#' posterior_interval(example_jm, regex_pars = c("Event", "Assoc"))
#' posterior_interval(example_jm, regexpars = "b[", prob = 0.95)
#' posterior_interval(example_jm, pars = "Long1|(Intercept)", prob = 0.5)
#' 
posterior_interval <- function(object, ...) {
  validate_stanjm_object(object)
  rstanarm:::posterior_interval(object = object, ...)
}
