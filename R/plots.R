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


#' Plot method for stanjm objects
#'
#' For models fit using \code{\link{stan_jm}}, there are a variety of
#' diagnostic and inference related plots that can be generated. The
#' user can specify which type of plot to display via the
#' \code{plotfun} argument.
#'
#' @method plot stanjm
#' @export
#' @templateVar stanjmArg x
#' @template args-stanjm-object
#' @template args-regex-pars
#' @param pars An optional character vector specifying a subset of parameters to
#'   display. Parameters can be specified by name or several shortcuts can be
#'   used.
#'   Using \code{pars = "long"} will display the
#'   parameter estimates for the longitudinal submodels only (excluding random
#'   effects, but including dispersion parameters).
#'   Using \code{pars = "event"} will display the
#'   parameter estimates for the event submodel only, including any association
#'   parameters.
#'   Using \code{pars = "assoc"} will display only the
#'   association parameters.
#'   Using \code{pars = "fixef"} will display all fixed effects, but not
#'   the random effects or the additional parameters such as dispersion, etc.
#'   Using \code{pars = "beta"} will display only the
#'   fixed effect regression coefficients (excluding the intercept terms).
#'   Using \code{pars = "alpha"} will display only the
#'   fixed effect intercept terms.
#'   The estimates for the random effects can be selected using \code{pars = "b"}
#'   or \code{pars = "varying"}.
#'   If both \code{pars} and \code{regex_pars} are set to \code{NULL} then all
#'   parameters are selected, including the log posterior.
#'   If both \code{pars} and \code{regex_pars} are set to \code{NULL} then all
#'   fixed effect regression coefficients are selected (with the exception
#'   of the intercept term if a Weibull baseline hazard was estimated).
#' @param plotfun A character string naming the plotting function to apply to
#'   the stanjm object. See \code{\link[rstanarm]{rstanarm-plots}} for the names and
#'   descriptions. Also see the \strong{Examples} section below. \code{plotfun} can be
#'   either the full name of the plotting function (e.g. \code{"stan_hist"}) or
#'   can be abbreviated to the part of the name following the underscore (e.g.
#'   \code{"hist"}). The default plot displays intervals and point estimates for
#'   the parameters.
#' @param ... Additional arguments to pass to \code{plotfun} (see
#'   \code{\link[rstanarm]{rstanarm-plots}}).
#'
#' @return A ggplot object (or several) that can be further
#'   customized using the \pkg{ggplot2} package.
#'
#' @seealso \code{\link[rstanarm]{rstanarm-plots}} for details on the individual plotting
#'   functions.
#'
#' @examples
#'
#' @importFrom rstan stan_plot stan_trace stan_scat stan_hist stan_dens stan_ac
#'   stan_diag stan_rhat stan_ess stan_mcse stan_par quietgg
#'
plot.stanjm <- function(x, plotfun = NULL, pars = NULL,
                         regex_pars = NULL, ...) {
  args <- set_plotting_args(x, pars, regex_pars, ...)
  fun <- set_plotting_fun(x, plotfun)
  do.call(fun, args)
}

# Prepare argument list to pass to plotting function
# @param x stanjm object
# @param pars,regex_pars user specified pars and regex_pars arguments (can be
#   missing)
set_plotting_args <- function(x, pars = NULL, regex_pars = NULL, ...) {
  args <- list(x, ...)
  M <- x$n_markers
  pars <- collect_pars(x, pars, regex_pars)

  # Collect parameter names
  nms <- collect_nms(rownames(x$stan_summary), M, value = TRUE)
  if (!is.null(pars)) {  # user specified pars
    pars2 <- NA
    if ("alpha" %in% pars) pars2 <- c(pars2, nms$alpha)
    if ("beta" %in% pars) pars2 <- c(pars2, nms$beta)
    if ("long" %in% pars) pars2 <- c(pars2, unlist(nms$y), unlist(nms$y_extra))
    if ("event" %in% pars) pars2 <- c(pars2, nms$e, nms$a, nms$e_extra)
    if ("assoc" %in% pars) pars2 <- c(pars2, nms$a)
    if ("fixef" %in% pars) pars2 <- c(pars2, unlist(nms$y), nms$e, nms$a)
    if ("b" %in% pars) pars2 <- c(pars2, nms$b)
    pars2 <- c(pars2, setdiff(pars,
                              c("alpha", "beta", "varying", "b",
                                "long", "event", "assoc", "fixef")))
    pars <- pars2[!is.na(pars2)]
  } else {  # no pars specified
    # Note: pars cannot be left equal to NULL since the
    # rstan plot functions will try to extract parameter names
    # using a call to
    #   names(object$coefficients)
    # which will throw up an error for 'stanjm' objects since it returns
    #   c("Long1", "Event")
    # and NOT the parameter names which are at the lower level of
    # the list!
    pars <- c(unlist(nms$y), nms$e, nms$a)
    pars <- setdiff(pars, "Event|(Intercept)")
  }

  args$pars <- check_plotting_pars(x, pars)
  return(args)
}

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
    STOP_sampling_only("pairs")
  requireNamespace("rstan")
  requireNamespace("KernSmooth")

  b <- b_names(rownames(x$stan_summary), value = TRUE)
  pars <- setdiff(rownames(x$stan_summary), b)
  dots <- list(...)
  user_pars <- dots[["pars"]]
  dots[["pars"]] <- NULL

  if (is.null(user_pars)) {  # user didn't specify pars
    pars <- drop_nms_with_spaces(pars)
    return(pairs(x$stanfit, pars = pars, ...))
  } else {  # user did specify pars
    if (any(user_pars %in% b))
      stop("pairs.stanjm does not yet allow group-level parameters in 'pars'.")
    user_pars <- drop_nms_with_spaces(user_pars)
    return(do.call("pairs", c(list(x$stanfit, pars = user_pars), dots)))
  }
}

# Function to remove variables from pars, if they have
# spaces in the variable name (since 'pairs' method seems to
# throw up an error if they are included)
# @param nms A character vector containing the names of variables
#   to be used in the 'pars' argument of the pairs.stanfit call
drop_nms_with_spaces <- function(nms) {
  sel <- grep(" ", nms)
  # pairs cannot handle variables with spaces in the variable
  # names -- the problem appears to be with passing these
  # variables in the 'pars' argument to the 'extract' function
  # (which is called within the 'pairs.stanfit' function)
  if (length(sel)) {
    cat("'pairs' cannot handle variable names with spaces. The ",
        "following variables will not be plotted: ",
        paste(nms[sel], collapse = ", "))
    nms <- nms[-sel]
  }
  nms
}
