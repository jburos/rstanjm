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

#' Plot estimated subject-specific longitudinal trajectory
#' 
#' Various plots comparing the observed outcome variable from one
#' of the longitudinal submodels \eqn{y} to simulated 
#' datasets \eqn{y^{rep}}{yrep} from the posterior predictive distribution.
#' This function is modelled on the \code{\link[rstanarm]{pp_check}} function
#' from the \pkg{rstanarm} package.
#' 
#' @method plot predict.stanjm
#' @export
#' 
#' @importFrom ggplot2 ggplot aes aes_string geom_line geom_smooth geom_ribbon 
#'   geom_point facet_wrap geom_vline labs ggplot_build theme_bw
plot.predict.stanjm <- function(object, ids = NULL, limits = c("ci", "pi", "none"), 
                                xlab = NULL, ylab = NULL, 
                                abline = FALSE, plot_observed = FALSE, facet_scales = "free_x", 
                                ci_geom_args = NULL, ...) {
  
  limits <- match.arg(limits)
  if (!(limits == "none")) ci <- (limits == "ci")
  y_var <- attr(object, "y_var")
  id_var <- attr(object, "id_var")
  time_var <- attr(object, "time_var")
  obs_dat <- attr(object, "observed_data")
  if (is.null(ylab)) ylab <- paste0("Long. response (", y_var, ")")
  if (is.null(xlab)) xlab <- paste0("Time (", time_var, ")")
  if (!is.null(ids)) {
    last_time <- attr(object, "last_time")[as.character(ids)]
    plot_dat <- object[object[[id_var]] %in% ids, , drop = FALSE]
    obs_dat <- obs_dat[obs_dat[[id_var]] %in% ids, , drop = FALSE]
  } else {
    ids <- attr(object, "ids")
    last_time <- attr(object, "last_time")
    plot_dat <- object
  }
  
  plot_dat$time <- plot_dat[[time_var]]
  plot_dat$id <- plot_dat[[id_var]]
  
  obs_dat$y <- obs_dat[[y_var]]
  obs_dat$time <- obs_dat[[time_var]]
  obs_dat$id <- obs_dat[[id_var]]
  
  geom_defaults <- list(color = "black", method = "loess", se = FALSE)
  geom_args <- rstanarm:::set_geom_args(geom_defaults, ...)

  lim_defaults <- list(alpha = 0.3)
  lim_args <- rstanarm:::set_geom_args(lim_defaults, ci_geom_args)

  obs_defaults <- list()
  obs_args <- rstanarm:::set_geom_args(obs_defaults)
  
  
  if (length(ids) > 60L) {
    stop("Too many individuals to plot for. Limit using 'ids' argument, or ",
         "estimate marginal longitudinal trajectories or those based on ",
         "specifying the 'newdata' argument.")
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
  if (abline) {
    graph_abline <- geom_vline(aes_string(xintercept = "last_time"), 
                               data.frame(id = ids, last_time = last_time), 
                               linetype = 2)
  } else graph_abline <- NULL
    
  
  graph + graph_limits + graph_obs + graph_abline + labs(x = xlab, y = ylab) 
}

