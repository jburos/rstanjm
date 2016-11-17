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

#' Plot the estimated subject-specific or marginal survival function
#' 
#' This generic \code{plot} method for \code{survfit.stanjm} objects will
#' plot the estimated subject-specific or marginal survival function
#' using the data frame returned by a call to \code{\link{posterior_survfit}}.
#' The call to \code{posterior_survfit} should ideally have included an
#' "extrapolation" of the survival function, obtained by setting the 
#' \code{extrapolate} argument to \code{TRUE}.
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
#' @param x A data frame and object of class \code{survfit.stanjm}
#'   returned by a call to the function \code{\link{posterior_survfit}}.
#'   The object contains point estimates and uncertainty interval limits
#'   for estimated values of the survival function.
#' @param limits A quoted character string specifying the type of limits to
#'   include in the plot. Can be one of: \code{"ci"} for the Bayesian
#'   posterior uncertainty interval for the estimated survival probability
#'   (often known as a credible interval); or \code{"none"} for no interval 
#'   limits.
#' @param ... Optional arguments passed to 
#'   \code{\link[ggplot2]{geom_line}} and used to control features
#'   of the plotted survival function.
#'      
#' @return A \code{ggplot} object, also of class \code{plot.survfit.stanjm}.
#'   This object can be further customised using the \pkg{ggplot2} package.
#'   It can also be passed to the function \code{\link{plot_stack}}.
#'   
#' @seealso \code{\link{posterior_survfit}}, \code{\link{plot_stack}},
#'   \code{\link{posterior_predict}}, \code{\link{plot.predict.stanjm}}      
#'   
#' @examples 
#' 
#'   # Run example model if not already loaded
#'   if (!exists("examplejm")) example(examplejm)
#'   
#'   # Obtain subject-specific conditional survival probabilities
#'   # for all individuals in the estimation dataset.
#'   ps1 <- posterior_survfit(examplejm, extrapolate = TRUE)
#'   
#'   # We then plot the conditional survival probabilities for
#'   # a subset of individuals
#'   plot(ps1, ids = c(7,13,16))
#'   
#'   # We can change or add attributes to the plot
#'   plot(ps1, ids = c(7,13,16), limits = "none")
#'   plot(ps1, ids = c(7,13,16), xlab = "Follow up time")
#'   plot(ps1, ids = c(7,13,16), ci_geom_args = list(fill = "red"),
#'        color = "blue", linetype = 2)
#'   plot(ps1, ids = c(7,13,16), facet_scales = "fixed")
#'   
#'   # Since the returned plot is also a ggplot object, we can
#'   # modify some of its attributes after it has been returned
#'   plot1 <- plot(ps1, ids = c(7,13,16))
#'   plot1 + 
#'     ggplot2::theme(strip.background = ggplot2::element_blank()) +
#'     ggplot2::coord_cartesian(xlim = c(0, 15)) +
#'     ggplot2::labs(title = "Some plotted survival functions")
#'     
#'   # We can also combine the plot(s) of the estimated 
#'   # subject-specific survival functions, with plot(s) 
#'   # of the estimated longitudinal trajectories for the
#'   # same individuals
#'   ps1 <- posterior_survfit(examplejm, ids = c(7,13,16))
#'   pt1 <- posterior_predict(examplejm, , ids = c(7,13,16),
#'                            interpolate = TRUE, extrapolate = TRUE)
#'   plot_surv <- plot(ps1) 
#'   plot_traj <- plot(pt1, vline = TRUE, plot_observed = TRUE)
#'   plot_stack(plot_traj, plot_surv)
#'    
#'   # Lastly, let us plot the standardised survival function
#'   # based on all individuals in our estimation dataset
#'   ps2 <- poserior_survfit(examplejm, standardise = TRUE, times = 0,
#'                           control = list(ext_points = 20))
#'   plot(ps2)   
#' 
#'    
plot.survfit.stanjm <- function(x, ids = NULL, 
                                limits = c("ci", "none"),  
                                xlab = NULL, ylab = NULL, facet_scales = "free", 
                                ci_geom_args = NULL, ...) {

  limits <- match.arg(limits)
  ci <- (limits == "ci")
  standardise <- attr(x, "standardise")
  id_var <- attr(x, "id_var")
  time_var <- attr(x, "time_var")
  if (is.null(xlab)) xlab <- paste0("Time (", time_var, ")")
  if (is.null(ylab)) ylab <- "Event free probability"
  if (!is.null(ids)) {
    if (standardise) 
      stop("'ids' argument cannot be specified when plotting standardised ",
           "survival probabilities.")
    if (!id_var %in% colnames(x))
      stop("Bug found: could not find 'id_var' column in the data frame.")
    ids_missing <- which(!ids %in% x[[id_var]])
    if (length(ids_missing))
      stop("The following 'ids' are not present in the survfit.stanjm object: ",
           paste(ids[[ids_missing]], collapse = ", "), call. = FALSE)
    x <- x[(x[[id_var]] %in% ids), , drop = FALSE]
  } else {
    ids <- if (!standardise) attr(x, "ids") else NULL
  }
  if (!standardise) x$id <- x[[id_var]]
  x$time <- x[[time_var]]
  
  geom_defaults <- list(color = "black")
  geom_args <- set_geom_args(geom_defaults, ...)  
  
  lim_defaults <- list(alpha = 0.3)
  lim_args <- do.call("set_geom_args", c(defaults = list(lim_defaults), ci_geom_args))
  
  if ((!standardise) && (length(ids) > 60L)) {
    stop("Too many individuals to plot for. Perhaps consider limiting ",
         "the number of individuals by specifying the 'ids' argument.")
  } else if ((!standardise) && (length(ids) > 1L)) {
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



