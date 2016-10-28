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
#
#' Graphical posterior predictive checks for the longitudinal submodel
#' 
#' Various plots comparing the observed outcome variable from one
#' of the longitudinal submodels \eqn{y} to simulated 
#' datasets \eqn{y^{rep}}{yrep} from the posterior predictive distribution.
#' This function is modelled on the \code{\link[rstanarm]{pp_check}} function
#' from the \pkg{rstanarm} package.
#' 
#' @export
#' @templateVar bdaRef (Ch. 6)
#' @templateVar stanjmArg object
#' @templateVar mArg m
#' @template reference-bda
#' @template args-stanjm-object
#' @template args-m
#' 
#' @param check The type of plot (possibly abbreviated) to show. One of 
#'   \code{"distributions"}, \code{"residuals"}, \code{"scatter"}, 
#'   \code{"test"}. See Details for descriptions.
#' @param nreps The number of \eqn{y^{rep}}{yrep} datasets to generate from the
#'   posterior predictive distribution (\code{\link{posterior_predict}}) and
#'   show in the plots. The default is \code{nreps=3} for
#'   \code{check="residuals"} and \code{nreps=8} for
#'   \code{check="distributions"}. If \code{check="test"}, \code{nreps} is
#'   ignored and the number of simulated datasets is the number of post-warmup
#'   draws from the posterior distribution. If \code{check="scatter"},
#'   \code{nreps} is not ignored but defaults to the number of post-warmup
#'   draws.
#' @param seed An optional \code{\link[=set.seed]{seed}} to pass to 
#'   \code{\link{posterior_predict}}.
#' @param overlay For \code{check="distributions"} only, should distributions be
#'   plotted as density estimates overlaid in a single plot (\code{TRUE}, the 
#'   default) or as separate histograms (\code{FALSE})?
#' @param test For \code{check="test"} only, a character vector (of length 1 or 
#'   2) naming a single function or a pair of functions. The function(s) should 
#'   take a vector input and return a scalar test statistic. See Details and
#'   Examples.
#' @param ... Optional arguments to geoms to control features of the plots 
#'   (e.g. \code{binwidth} if the plot is a histogram).
#' 
#' @return A ggplot object that can be further customized using the
#'   \pkg{ggplot2} package.
#'   
#' @details 
#' Note that the posterior predictive checks here condition on all 
#' observed data, that is for a predicted observation \eqn{y^{rep}_{i}(t)}
#' the predictive distribution contains information from data that may
#' have been observed later in time, that is some time \eqn{s} where 
#' \eqn{s>t}. The predictions are therefore suitable for assessing model
#' fit, but should not be used to dynamic predictions where we want to 
#' condition only on previous measurements for the individual.\cr
#' \cr
#' Descriptions of the plots corresponding to the different values of 
#' \code{check}:
#' \describe{
#'  \item{\code{distributions}}{The distributions of \eqn{y} and \code{nreps} 
#'  simulated \eqn{y^{rep}}{yrep} datasets.} 
#'  \item{\code{residuals}}{The distributions of residuals computed from \eqn{y}
#'  and each of \code{nreps} simulated datasets. For binomial data, binned 
#'  residual plots are generated (similar to \code{\link[arm]{binnedplot}}).}
#'  \item{\code{scatter}}{If \code{nreps} is \code{NULL} then \eqn{y} is plotted
#'  against the average values of \eqn{y^{rep}}{yrep}, i.e., the points
#'  \eqn{(y_n, \bar{y}^{rep}_n),\, n = 1, \dots, N}{(y_n, mean(yrep_n)), n =
#'  1,...,N}, where each \eqn{y^{rep}_n}{yrep_n} is a vector of length equal to
#'  the number of posterior draws. If \code{nreps} is a (preferably small)
#'  integer, then only \code{nreps} \eqn{y^{rep}}{yrep} datasets are simulated
#'  and they are each plotted separately against \eqn{y}.}
#'  \item{\code{test}}{The distribution of a single test statistic
#'  \eqn{{T(y^{rep})}}{T(yrep)} or a pair of test statistics over the
#'  \code{nreps} simulated datasets. If the \code{test} argument specifies only
#'  one function then the resulting plot is a histogram of
#'  \eqn{{T(y^{rep})}}{T(yrep)} and the value of the test statistic in the 
#'  observed data, \eqn{T(y)}, is shown in the plot as a vertical line. If two 
#'  functions are specified then the plot is a scatterplot and \eqn{T(y)} is 
#'  shown as a large point.}
#' }
#' 
#' @note For binomial data, plots of \eqn{y} and \eqn{y^{rep}}{yrep} show the
#'   proportion of 'successes' rather than the raw count.
#' 
#' @seealso \code{\link{posterior_predict}} for drawing from the posterior 
#'   predictive distribution for the longitudinal submodel, 
#'   \code{\link{posterior_survfit}} for the estimated survival function based
#'   on the posterior distribution of the model parameters, and 
#'   \code{\link{ps_check}} for graphical checks of the estimated survival 
#'   function.
#' 
#' @examples 
#' 
pp_check <- function(object, m = 1, check = "distributions", nreps = NULL, 
                     seed = NULL, overlay = TRUE, test = "mean", ...) {
  validate_stanjm_object(object)
  if (rstanarm:::used.optimizing(object)) 
    rstanarm:::STOP_not_optimizing("pp_check")
  
  checks <- c("distributions", "residuals", "scatter", "test")
  fn <- switch(match.arg(arg = check, choices = checks),
               'distributions' = "pp_check_dist",
               'residuals' = "pp_check_resid",
               'test' = "pp_check_stat",
               'scatter' = "pp_check_scatter")
  if (is.null(nreps) && !fn %in%  c("pp_check_stat", "pp_check_scatter"))
    nreps <- ifelse(fn == "pp_check_dist", 8, 3)
  if (!is.null(nreps) && fn == "pp_check_stat") {
    warning("'nreps' is ignored if check='test'", call. = FALSE)
    nreps <- NULL
  }

  is_binomial <- if (is(object, "polr") && !rstanarm:::is_scobit(object))
    FALSE else rstanarm:::is.binomial(family(object)[[m]]$family)
  if (is_binomial && fn == "pp_check_resid") {
    graph <- rstanarm:::pp_check_binned_resid(object, n = nreps, ...)
    return(graph)
  }
  
  y <- get_y(object)[[m]]
  yrep <- posterior_predict(object, m, draws = nreps, seed = seed)
  if (is(object, "polr")) {
    y <- as.integer(y)
    yrep <- apply(yrep, 2, function(x) as.integer(as.factor(x)))
  }
  if (is_binomial) {
    if (NCOL(y) == 2L) {
      trials <- rowSums(y)
      y <- y[, 1L] / trials
      yrep <- sweep(yrep, 2L, trials, "/")
    } else if (is.factor(y))
      y <- rstanarm:::fac2bin(y)
  }
  
  if (fn == "pp_check_dist") {
    args <- list(y = y, yrep = yrep, n = nreps, overlay = overlay, ...)
  } else if (fn == "pp_check_resid") {
    args <- list(y = y, yrep = yrep, n = nreps, ...)
  } else if (fn == "pp_check_stat") {
    args <- list(y = y, yrep = yrep, test = test, ...)
  } else if (fn == "pp_check_scatter") {
    args <- list(y = y, yrep = yrep, n = nreps, ...)
  }
  
  fn <- utils::getFromNamespace(fn, "rstanarm")
  do.call(fn, args)
}
