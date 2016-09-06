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
plot.survfit.stanjm <- function(object, ci = TRUE, limits = c(.025, .975), ids = NULL, ...) {

  marginalised <- attr(object, "marginalised")
  n_increments <- attr(object, "n_increments")
  if (!is.null(ids)) {
    times <- lapply(object$times, function(x) 
      x[names(x) %in% ids])
    survprobs <- lapply(object$survprobs, function(x) 
      x[, colnames(x) %in% ids, drop = FALSE])
  } else {
    ids <- attr(object, "ids")
    times <- object$times
    survprobs <- object$survprobs
  }
  
  times <- do.call(cbind, times)
  if (NROW(times) == 1L) 
    times <- as.vector(times)  
  survprob_median <- sapply(survprobs, function(x) apply(x, 2, median))
  survprob_lb <- sapply(survprobs, function(x) apply(x, 2, quantile, limits[1])) 
  survprob_ub <- sapply(survprobs, function(x) apply(x, 2, quantile, limits[2])) 
  
  defaults <- list(color = "black")
  geom_args <- rstanarm:::set_geom_args(defaults, ...)  
  
  if ((!marginalised) && (length(ids) > 60L)) {
    stop("Too many individuals to plot for. Limit using 'ids' argument, or ",
         "estimate marginal survival probabilities by using 'posterior_survfit' ",
         "with 'marginalised' set to TRUE.")
  } else if ((!marginalised) && (length(ids) > 1L)) {
    dats <- list(times, survprob_median, survprob_lb, survprob_ub)
    v_names <- c("times", "med", "lb", "ub")
    plot_dats <- mapply(prep_for_facetwrap, dats,
                        v_names, n_increments, SIMPLIFY = FALSE)
    plot_dat <- Reduce(function(...) merge(...), plot_dats)
    graph <- ggplot(plot_dat, aes_string(x = "times", y = "med")) +
      theme_bw() +
      do.call("geom_line", geom_args) +
      coord_cartesian(ylim = c(0, 1)) +      
      labs(x = "Time", y = "Event free probability") +
      facet_wrap(~ id, scales = "free")
    if (ci) {
      graph_limits <- geom_ribbon(aes_string(ymin = "lb", ymax = "ub"), 
                                  alpha = 0.3)      
    } else graph_limits <- NULL
  } else {
    plot_dat <- data.frame(times = times, med = survprob_median,
                           lb = survprob_lb, ub = survprob_ub)
    graph <- ggplot(plot_dat, aes_string(x = "times", y = "med")) + 
      theme_bw() +
      do.call("geom_line", geom_args) + 
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Time", y = "Event free probability")
    if (ci) {
      graph_limits <- geom_ribbon(aes_string(ymin = "lb", ymax = "ub"), 
                                  alpha = 0.3)      
    } else graph_limits <- NULL

  }    

  graph + graph_limits
}


# reshape data for facet_wrap
prep_for_facetwrap <- function(x, v_names, n_increments, ...) {
  x <- as.data.frame(x)
  x$id <- rownames(x)
  x <- reshape(x, direction = "long", varying = paste0("increment", 0:n_increments), 
               v.names = v_names, timevar = "obs", idvar = "id", ...)
  x
}


#-------------------------------------------------------------

# default plotting attributes
.PP_FILL <- "skyblue"
.PP_DARK <- "skyblue4"
.PP_VLINE_CLR <- "#222222"
.PP_YREP_CLR <- "#487575"
.PP_YREP_FILL <- "#222222"



