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

#' Plot the estimated subject-specific or marginal longitudinal trajectory
#' 
#' This generic \code{plot} method for \code{predict.stanjm} objects will
#' plot the estimated subject-specific or marginal longitudinal trajectory
#' using the data frame returned by a call to \code{\link{posterior_predict}}.
#' To ensure that enough data points are available to plot the longitudinal
#' trajectory, it is assumed that the call to \code{\link{posterior_predict}}
#' would have included the argument \code{interpolate = TRUE}, and perhaps
#' also \code{extrapolate = TRUE} (the latter being optional, depending on 
#' whether or not the user wants to see extrapolation of the longitudinal 
#' trajectory beyond the last observation time).
#' 
#' @method plot predict.stanjm
#' @export
#' @importFrom ggplot2 ggplot aes aes_string geom_line geom_smooth geom_ribbon 
#'   geom_point facet_wrap geom_vline labs ggplot_build theme_bw
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
#' @param x A data frame and object of class \code{predict.stanjm}
#'   returned by a call to the function \code{\link{posterior_predict}}.
#'   The object contains point estimates and uncertainty interval limits
#'   for the fitted values of the longitudinal response.
#' @param limits A quoted character string specifying the type of limits to
#'   include in the plot. Can be one of: \code{"ci"} for the Bayesian
#'   posterior uncertainty interval for the estimated mean longitudinal
#'   response (often known as a credible interval);
#'   \code{"pi"} for the prediction interval for the estimated (raw)
#'   longitudinal response; or \code{"none"} for no interval limits.
#' @param vline A logical. If \code{TRUE} then a vertical dashed line
#'   is added to the plot indicating the event or censoring time for
#'   the individual. Can only be used if each plot within the figure
#'   is for a single individual.
#' @param plot_observed A logical. If \code{TRUE} then the observed
#'   longitudinal measurements are overlaid on the plot.
#' @param ... Optional arguments passed to 
#'   \code{\link[ggplot2]{geom_smooth}} and used to control features
#'   of the plotted longitudinal trajectory.
#'   
#' @return A \code{ggplot} object, also of class \code{plot.predict.stanjm}.
#'   This object can be further customised using the \pkg{ggplot2} package.
#'   It can also be passed to the function \code{\link{plot_stack}}.
#'   
#' @seealso \code{\link{posterior_predict}}, \code{\link{plot_stack}},
#'   \code{\link{posterior_survfit}}, \code{\link{plot.survfit.stanjm}}   
#'     
#' @examples 
#' \donttest{
#' }
#' 
plot.predict.stanjm <- function(x, ids = NULL, limits = c("ci", "pi", "none"), 
                                xlab = NULL, ylab = NULL, vline = FALSE, 
                                plot_observed = FALSE, facet_scales = "free_x", 
                                ci_geom_args = NULL, ...) {
  
  limits <- match.arg(limits)
  if (!(limits == "none")) ci <- (limits == "ci")
  y_var <- attr(x, "y_var")
  id_var <- attr(x, "id_var")
  time_var <- attr(x, "time_var")
  obs_dat <- attr(x, "observed_data")
  if (is.null(ylab)) ylab <- paste0("Long. response (", y_var, ")")
  if (is.null(xlab)) xlab <- paste0("Time (", time_var, ")")
  if (!id_var %in% colnames(x))
    stop("Bug found: could not find 'id_var' column in the data frame.")
  if (!is.null(ids)) {
    ids_missing <- which(!ids %in% x[[id_var]])
    if (length(ids_missing))
      stop("The following 'ids' are not present in the predict.stanjm object: ",
           paste(ids[[ids_missing]], collapse = ", "), call. = FALSE)
    plot_dat <- x[x[[id_var]] %in% ids, , drop = FALSE]
    obs_dat <- obs_dat[obs_dat[[id_var]] %in% ids, , drop = FALSE]
  } else {
    plot_dat <- x
  }
  
  # Reorder ids to match order in plotting data
  ids <- factor(unique(plot_dat[[id_var]]))
  last_time <- attr(x, "last_time")[as.character(ids)]
  
  plot_dat$time <- plot_dat[[time_var]]
  plot_dat$id <- plot_dat[[id_var]]
  
  obs_dat$y <- obs_dat[[y_var]]
  obs_dat$time <- obs_dat[[time_var]]
  obs_dat$id <- obs_dat[[id_var]]
  
  geom_defaults <- list(color = "black", method = "loess", se = FALSE)
  geom_args <- set_geom_args(geom_defaults, ...)

  lim_defaults <- list(alpha = 0.3)
  lim_args <- set_geom_args(lim_defaults, ci_geom_args)

  obs_defaults <- list()
  obs_args <- set_geom_args(obs_defaults)
  
  
  if (length(ids) > 60L) {
    stop("Too many individuals to plot for. Perhaps limit the number of ",
         "individuals by specifying the 'ids' argument.")
  } else if (length(ids) > 1L) {
    geom_mapp <- list(mapping = aes_string(x = "time", y = "yfit"), 
                      data = plot_dat)
    graph <- ggplot() + theme_bw() +
               do.call("geom_smooth", c(geom_mapp, geom_args)) +
               facet_wrap(~ id, scales = facet_scales)
    if (!(limits == "none")) {
      graph_smoothlim <- ggplot(plot_dat) + 
        geom_smooth(aes_string(x = "time", y = if (ci) "ci_lb" else "pi_lb"), 
                    method = "loess", se = FALSE) +
        geom_smooth(aes_string(x = "time", y = if (ci) "ci_ub" else "pi_ub"), 
                    method = "loess", se = FALSE) +
        facet_wrap(~ id, scales = facet_scales)
      build_smoothlim <- ggplot_build(graph_smoothlim)
      df_smoothlim <- data.frame(PANEL = build_smoothlim$data[[1]]$PANEL,
                                 time = build_smoothlim$data[[1]]$x,
                                 lb = build_smoothlim$data[[1]]$y,
                                 ub = build_smoothlim$data[[2]]$y)
      panel_id_map <- build_smoothlim$panel$layout[, c("PANEL", "id"), drop = FALSE]
      df_smoothlim <- merge(df_smoothlim, panel_id_map)
      lim_mapp <- list(mapping = aes_string(x = "time", ymin = "lb", ymax = "ub"), 
                       data = df_smoothlim)
      graph_limits <- do.call("geom_ribbon", c(lim_mapp, lim_args))
    } else graph_limits <- NULL
  } else {
    geom_mapp <- list(mapping = aes_string(x = "time", y = "yfit"), 
                      data = plot_dat)
    graph <- ggplot() + theme_bw() + 
               do.call("geom_smooth", c(geom_mapp, geom_args))
    if (!(limits == "none")) {
      graph_smoothlim <- ggplot(plot_dat) + 
        geom_smooth(aes_string(x = "time", y = if (ci) "ci_lb" else "pi_lb"), 
                    method = "loess", se = FALSE) +
        geom_smooth(aes_string(x = "time", y = if (ci) "ci_ub" else "pi_ub"), 
                    method = "loess", se = FALSE)
      build_smoothlim <- ggplot_build(graph_smoothlim)
      df_smoothlim <- data.frame(time = build_smoothlim$data[[1]]$x,
                                 lb = build_smoothlim$data[[1]]$y,
                                 ub = build_smoothlim$data[[2]]$y) 
      lim_mapp <- list(mapping = aes_string(x = "time", ymin = "lb", ymax = "ub"), 
                       data = df_smoothlim)
      graph_limits <- do.call("geom_ribbon", c(lim_mapp, lim_args))
    } else graph_limits <- NULL
  }    
  if (plot_observed) {
    obs_mapp <- list(mapping = aes_string(x = "time", y = "y"), 
                     data = obs_dat)
    graph_obs <- do.call("geom_point", c(obs_mapp, obs_args)) 
  } else graph_obs <- NULL
  if (vline) {
    graph_vline <- geom_vline(aes_string(xintercept = "last_time"), 
                               data.frame(id = ids, last_time = last_time), 
                               linetype = 2)
  } else graph_vline <- NULL
    
  
  ret <- graph + graph_limits + graph_obs + graph_vline + labs(x = xlab, y = ylab) 
  class_ret <- class(ret)
  class(ret) <- c("plot.predict.stanjm", class_ret)
  ret
  
}

