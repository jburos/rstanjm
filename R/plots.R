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
#' @templateVar stanjmArg x
#' @template args-stanjm-object
#' @param ... Arguments to pass to \code{\link[rstan]{pairs.stanfit}}. 
#' 
#' @description See \code{\link[rstan]{pairs.stanfit}} for details.
#' @details See the Details section in \code{\link[rstan]{pairs.stanfit}}.
#' @importFrom graphics pairs
#' @examples
#' if (!exists("examplejm")) example(examplejm)
#' pairs(examplejm, pars = c("Long1|(Intercept)", "log-posterior"))
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

