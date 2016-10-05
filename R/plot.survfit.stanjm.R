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
#' 
#' @importFrom ggplot2 ggplot aes_string geom_line geom_ribbon 
#'   facet_wrap labs coord_cartesian
plot.survfit.stanjm <- function(object, ids = NULL, ci = TRUE,  
                                xlab = NULL, ylab = NULL, ...) {

  marginalised <- attr(object, "marginalised")
  n_increments <- attr(object, "n_increments")
  id_var <- attr(object, "id_var")
  time_var <- attr(object, "time_var")
  if (is.null(xlab)) xlab <- paste0("Time (", time_var, ")")
  if (is.null(ylab)) ylab <- "Event free probability"
  if (!is.null(ids)) {
    ids_missing <- which(!ids %in% object[[id_var]])
    if (length(ids_missing))
      stop("The following 'ids' are not present in the survfit.stanjm object: ",
           paste(ids[[ids_missing]], collapse = ", "), call. = FALSE)
    object <- object[(object[[id_var]] %in% ids), , drop = FALSE]
  } else {
    ids <- if (!marginalised) attr(object, "ids") else NULL
  }
  object$id <- object[[id_var]]
  object$time <- object[[time_var]]
  
  defaults <- list(color = "black")
  geom_args <- rstanarm:::set_geom_args(defaults, ...)  
  
  if ((!marginalised) && (length(ids) > 60L)) {
    stop("Too many individuals to plot for. Limit using 'ids' argument, or ",
         "estimate marginal survival probabilities by using 'posterior_survfit' ",
         "with 'marginalised' set to TRUE.")
  } else if ((!marginalised) && (length(ids) > 1L)) {
    graph <- ggplot(object, aes_string(x = "time", y = "survpred")) +
      theme_bw() +
      do.call("geom_line", geom_args) +
      coord_cartesian(ylim = c(0, 1)) +      
      facet_wrap(~ id, scales = "free")
    if (ci) {
      graph_limits <- geom_ribbon(aes_string(ymin = "ci_lb", ymax = "ci_ub"), 
                                  alpha = 0.3)      
    } else graph_limits <- NULL
  } else {
    graph <- ggplot(object, aes_string(x = "time", y = "survpred")) + 
      theme_bw() +
      do.call("geom_line", geom_args) + 
      coord_cartesian(ylim = c(0, 1))
    if (ci) {
      graph_limits <- geom_ribbon(aes_string(ymin = "ci_lb", ymax = "ci_ub"), 
                                  alpha = 0.3)      
    } else graph_limits <- NULL

  }    

  graph + graph_limits + labs(x = xlab, y = ylab) 
}


#-------------------------------------------------------------

# default plotting attributes
.PP_FILL <- "skyblue"
.PP_DARK <- "skyblue4"
.PP_VLINE_CLR <- "#222222"
.PP_YREP_CLR <- "#487575"
.PP_YREP_FILL <- "#222222"



