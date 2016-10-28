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

#' Plot estimated subject-specific or marginal survival function
#' 
#' Various plots comparing the observed outcome variable from one
#' of the longitudinal submodels \eqn{y} to simulated 
#' datasets \eqn{y^{rep}}{yrep} from the posterior predictive distribution.
#' This function is modelled on the \code{\link[rstanarm]{pp_check}} function
#' from the \pkg{rstanarm} package.
#'    
#' @method plot survfit.stanjm
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_line geom_ribbon 
#'   facet_wrap labs coord_cartesian
#'   
#' @templateVar idsArg ids
#' @templateVar labsArg xlab,ylab
#' @templateVar scalesArg facet_scales
#' @templateVar cigeomArg ci_geom_args
#' @template args-ids
#' @template args-labs
#' @template args-scales
#' @template args-ci-geom-args
#'  
#' @param x A data frame and x of class \code{survfit.stanjm}
#'   returned by a call to the function \code{\link{posterior_survfit}}.
#'   The x contains point estimates and uncertainty interval limits
#'   for estimated values of the survival function.
#' @param limits A quoted character string specifying the type of limits to
#'   include in the plot. Can be one of: \code{"ci"} for the Bayesian
#'   posterior uncertainty interval (often known as a credible interval);
#'   or \code{"none"} for no interval limits.
#' @param ... Optional arguments passed to 
#'   \code{\link[ggplot2]{geom_line}} and used to control features
#'   of the plotted trajectory.
#'      
#' @return A \code{ggplot} x, also of class \code{plot.survfit.stanjm}.
#'   This x can be further customised using the \pkg{ggplot2} package.
#'   It can also be passed to the function \code{\link{plot_stack}}.
#'   
#' @seealso \code{\link{posterior_survfit}}, \code{\link{plot_stack}}       
#'   
#' @examples 
#'    
plot.survfit.stanjm <- function(x, ids = NULL, 
                                limits = c("ci", "none"),  
                                xlab = NULL, ylab = NULL, facet_scales = "free", 
                                ci_geom_args = NULL, ...) {

  marginalised <- attr(x, "marginalised")
  n_increments <- attr(x, "n_increments")
  limits <- match.arg(limits)
  ci <- (limits == "ci")
  id_var <- attr(x, "id_var")
  time_var <- attr(x, "time_var")
  if (is.null(xlab)) xlab <- paste0("Time (", time_var, ")")
  if (is.null(ylab)) ylab <- "Event free probability"
  if (!is.null(ids)) {
    ids_missing <- which(!ids %in% x[[id_var]])
    if (length(ids_missing))
      stop("The following 'ids' are not present in the survfit.stanjm x: ",
           paste(ids[[ids_missing]], collapse = ", "), call. = FALSE)
    x <- x[(x[[id_var]] %in% ids), , drop = FALSE]
  } else {
    ids <- if (!marginalised) attr(x, "ids") else NULL
  }
  x$id <- x[[id_var]]
  x$time <- x[[time_var]]
  
  geom_defaults <- list(color = "black")
  geom_args <- rstanarm:::set_geom_args(geom_defaults, ...)  
  
  lim_defaults <- list(alpha = 0.3)
  lim_args <- rstanarm:::set_geom_args(lim_defaults, ci_geom_args)
  
  if ((!marginalised) && (length(ids) > 60L)) {
    stop("Too many individuals to plot for. Limit using 'ids' argument, or ",
         "estimate marginal survival probabilities by using 'posterior_survfit' ",
         "with 'marginalised' set to TRUE.")
  } else if ((!marginalised) && (length(ids) > 1L)) {
    graph <- ggplot(x, aes_string(x = "time", y = "survpred")) +
      theme_bw() +
      do.call("geom_line", geom_args) +
      coord_cartesian(ylim = c(0, 1)) +      
      facet_wrap(~ id, scales = facet_scales)
    if (ci) {
      lim_mapp <- list(mapping = aes_string(ymin = "ci_lb", ymax = "ci_ub"))
      graph_limits <- do.call("geom_ribbon", c(lim_mapp, lim_args))
    } else graph_limits <- NULL
  } else {
    graph <- ggplot(x, aes_string(x = "time", y = "survpred")) + 
      theme_bw() +
      do.call("geom_line", geom_args) + 
      coord_cartesian(ylim = c(0, 1))
    if (ci) {
      lim_mapp <- list(mapping = aes_string(ymin = "ci_lb", ymax = "ci_ub"))
      graph_limits <- do.call("geom_ribbon", c(lim_mapp, lim_args))
    } else graph_limits <- NULL

  }    

  ret <- graph + graph_limits + labs(x = xlab, y = ylab) 
  class_ret <- class(ret)
  class(ret) <- c("plot.survfit.stanjm", class_ret)
  ret
}


#-------------------------------------------------------------

# default plotting attributes
.PP_FILL <- "skyblue"
.PP_DARK <- "skyblue4"
.PP_VLINE_CLR <- "#222222"
.PP_YREP_CLR <- "#487575"
.PP_YREP_FILL <- "#222222"



