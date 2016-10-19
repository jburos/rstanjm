# Part of the rstanarm package for estimating model parameters
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


#' Pairs method for stanjm objects
#' 
#' @method pairs stanjm
#' @export
#' @param x A stanjm object.
#' @param ... Arguments to pass to \code{\link[rstan]{pairs.stanfit}}. 
#' 
#' @description See \code{\link[rstan]{pairs.stanfit}} for details.
#' @details See the Details section in \code{\link[rstan]{pairs.stanfit}}.
#' @importFrom graphics pairs
#' @examples
#' if (!exists("examplejm")) example(examplejm)
#' pairs(examplejm, pars = c("(Intercept)", "log-posterior"))
#' 
pairs.stanjm <- function(x, ...) {
  if (!rstanarm:::used.sampling(x)) 
    rstanarm:::STOP_sampling_only("pairs")
  requireNamespace("rstan")
  requireNamespace("KernSmooth")
  
  if (rstanarm:::is.mer(x)) {
    b <- b_names(rownames(x$stan_summary), value = TRUE)
    pars <- setdiff(rownames(x$stan_summary), b)
    pars <- setdiff(pars, c("Assoc|Long1:eta value", "Event|weibull shape"))
    dots <- list(...)
    if (is.null(dots[["pars"]]))
      return(pairs(x$stanfit, pars = pars, ...)) 
    
    if (any(dots[["pars"]] %in% b))
      stop("pairs.stanjm does not yet allow group-level parameters in 'pars'.")
  }
  
  pairs(x$stanfit, ...)
}

#' Plots for rstanarm models
#' 
#' Models fit using \code{algorithm='sampling'}, \code{"meanfield"}, or
#' \code{"fullrank"} are compatible with a variety of plotting functions from
#' the \pkg{rstan} package. Each function returns at least one
#' \code{\link[ggplot2]{ggplot}} object that can be customized further using the
#' \pkg{ggplot2} package. The plotting functions described here can be called
#' using the \code{\link[=plot.stanreg]{plot method}} for stanreg objects 
#' without loading the \pkg{rstan} package. For example usage see 
#' \code{\link{plot.stanreg}}.
#' 
#' @name rstanarm-plots
#' 
#' @section Plotting functions:
#' 
#' \describe{
#' \item{Posterior intervals and point
#' estimates}{\code{\link[rstan]{stan_plot}}}
#' \item{Traceplots}{\code{\link[rstan]{stan_trace}}}
#' \item{Histograms}{\code{\link[rstan]{stan_hist}}}
#' \item{Kernel density estimates}{\code{\link[rstan]{stan_dens}}}
#' \item{Scatterplots}{\code{\link[rstan]{stan_scat}}}
#' \item{Diagnostics for Hamiltonian Monte Carlo and the No-U-Turn
#' Sampler}{\code{\link[rstan]{stan_diag}}}
#' \item{Rhat}{\code{\link[rstan]{stan_rhat}}}
#' \item{Ratio of effective sample size to total posterior sample
#' size}{\code{\link[rstan]{stan_ess}}}
#' \item{Ratio of Monte Carlo standard error to posterior standard
#' deviation}{\code{\link[rstan]{stan_mcse}}}
#' \item{Autocorrelation}{\code{\link[rstan]{stan_ac}}}
#' }
#' 
#' @seealso \code{\link{plot.stanreg}} for how to call the \code{plot} method, 
#'   \code{\link{shinystan}} for interactive model exploration,
#'   \code{\link{pp_check}} for graphical posterior predicive checking.
#'
#' @examples
#' # See examples at help("plot.stanreg", package = "rstanarm")
NULL
