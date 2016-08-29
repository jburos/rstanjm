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

#' Graphical posterior predictive checks for the event submodel
#' 
#' Various plots comparing the observed outcome variable from one
#' of the longitudinal submodels \eqn{y} to simulated 
#' datasets \eqn{y^{rep}}{yrep} from the posterior predictive distribution.
#' This function is modelled on the \code{\link[rstanarm]{pp_check}} function
#' from the \pkg{rstanarm} package.
#' 
#' @export
#' 
#' @importFrom ggplot2 ggplot aes_string geom_line labs facet_wrap
ps_check <- function(object, type = "survival", ci_limit = c(.025, .975), 
                     draws = NULL, seed = NULL, ...) {
  validate_stanjm_object(object)

  dat <- posterior_survfit(object, marginalised = TRUE, 
                           draws = draws, seed = seed, ...)
  
  if (type == "survival") {
    
    times <- do.call(cbind, dat$times)
    if (NROW(times) == 1L) 
      times <- as.vector(times)
    survprob_median <- sapply(dat$survprobs, function(x) apply(x, 2, median)) 
    survprob_lb <- sapply(dat$survprobs, function(x) apply(x, 2, quantile, ci_limit[1])) 
    survprob_ub <- sapply(dat$survprobs, function(x) apply(x, 2, quantile, ci_limit[2])) 
 
    defaults <- list(color = "black")
    geom_args <- rstanarm:::set_geom_args(defaults, ...)  
       
    plot_dat <- data.frame(x = times, y = survprob_median)
    graph <- ggplot(plot_dat, aes_string("x", "y")) + 
      do.call("geom_line", geom_args) +
      labs(x = "Time", y = "Event free probability")
    
    # overlay KM curve
    
  }
  graph
}


# default plotting attributes
.PP_FILL <- "skyblue"
.PP_DARK <- "skyblue4"
.PP_VLINE_CLR <- "#222222"
.PP_YREP_CLR <- "#487575"
.PP_YREP_FILL <- "#222222"



