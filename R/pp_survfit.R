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
#' 
#' @importFrom ggplot2 ggplot aes_string geom_line labs facet_wrap
pp_survfit <- function(object, type = "survival", ids, draws = NULL,
                      ci_limit = c(.025, .975), 
                      newdataEvent = NULL, newdataLong = NULL, seed = NULL, ...) {
  validate_stanjm_object(object)

  dat <- posterior_predictEvent(object, newdataEvent, newdataLong, ids,
                                draws, seed, ...)
  
  if (type == "survival") {
    
    times <- do.call(cbind, dat$times)
    if (NROW(times) == 1L) 
      times <- as.vector(times)
    survprob_median <- sapply(dat$survprobs, function(x) apply(x, 2, median)) 
    survprob_lb <- sapply(dat$survprobs, function(x) apply(x, 2, quantile, ci_limit[1])) 
    survprob_ub <- sapply(dat$survprobs, function(x) apply(x, 2, quantile, ci_limit[2])) 
 
    defaults <- list(color = "black")
    geom_args <- rstanarm:::set_geom_args(defaults, ...)  
       
    if (length(ids) > 1){
      times <- prep_for_facetwrap(times, v.names = "times")
      survprob_median <- prep_for_facetwrap(survprob_median, v.names = "survprob_median")
      survprob_lb <- prep_for_facetwrap(survprob_lb, v.names = "survprob_lb")
      survprob_ub <- prep_for_facetwrap(survprob_ub, v.names = "survprob_ub")
      plot_dat <- Reduce(function(...) merge(...), 
                         list(times, survprob_median, survprob_lb, survprob_ub))
      plot_dat <<- plot_dat
      graph <- ggplot(plot_dat, aes_string("times", "survprob_median")) + 
        do.call("geom_line", geom_args) +
        labs(x = "Time", y = "Event free probability") +
        facet_wrap(~ id, scales = "free")
    } else {
      plot_dat <- data.frame(x = times, y = survprob_median)
      graph <- ggplot(plot_dat, aes_string("x", "y")) + 
        do.call("geom_line", geom_args) +
        labs(x = "Time", y = "Event free probability")
    }
        
  }
  graph
}


# default plotting attributes
.PP_FILL <- "skyblue"
.PP_DARK <- "skyblue4"
.PP_VLINE_CLR <- "#222222"
.PP_YREP_CLR <- "#487575"
.PP_YREP_FILL <- "#222222"


# reshape data for facet_wrap
prep_for_facetwrap <- function(x, ...) {
  x <- as.data.frame(x)
  x$id <- rownames(x)
  x <- reshape(x, direction = "long", varying = paste0("V", 1:(NCOL(x)-1)), 
               timevar = "obs", idvar = "id", ...)
  x
}


