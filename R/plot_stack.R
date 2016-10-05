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

#' 
#' @export
#' 
#' @importFrom ggplot2 ggplot_build facet_wrap aes_string expand_limits
#' 
plot_stack <- function(yplot, survplot) {
    
  if (!is(yplot, "list")) yplot <- list(yplot)
  
  y_build <- lapply(yplot, ggplot_build)
  y_layout <- lapply(y_build, function(x) x$panel$layout)
  y_ids <- lapply(y_layout, function(x)
    if (!"id" %in% colnames(x)) NULL else x[["id"]])
  
  e_build <- ggplot_build(survplot)
  e_layout <- e_build$panel$layout    
  e_ids <- if (!"id" %in% colnames(e_layout)) NULL else e_layout[["id"]]
  
  lapply(y_ids, function(x, e_ids) {
    if (!identical(x, e_ids)) 
      stop("The individuals in the 'yplot' and 'survplot' appear to differ. Please ",
           "reestimate the plots using a common 'ids' argument.", call. = FALSE)
    }, e_ids = e_ids)
  
  vline <- lapply(seq_along(y_build), function(m) {
    L <- length(y_build[[m]]$data)
    dat <- y_build[[m]]$data[[L]]
    if (!"xintercept" %in% colnames(dat)) {
      found <- FALSE
    } else {
      found <- TRUE
      dat <- dat[, c("PANEL", "xintercept"), drop = FALSE] 
      if (NROW(y_layout[[m]]) > 1) {
        panel_id_map <- y_layout[[m]][, c("PANEL", "id"), drop = FALSE]
        dat <- merge(dat, panel_id_map, by = "PANEL")
      }
      dat <- dat[, grep("PANEL", colnames(dat), invert = TRUE), drop = FALSE]
      colnames(dat) <- gsub("xintercept", paste0("xintercept", m), colnames(dat), fixed = TRUE)
    }
    list(dat = dat, found = found)
  })
  vline_found <- any(sapply(vline, function(x) x$found))
  if (!vline_found)
    warning("Could not find vertical line indicating last observation time in the ",
            "plots of the longitudinal trajectory. You may wish to plot the longitudinal ",
            "trajectories again with 'abline = TRUE'.", immediate. = TRUE)
  vline_dat <- lapply(vline, function(x) x$dat)
  vline_alldat <- Reduce(function(...) merge(..., all = TRUE), vline_dat)
  vline_alldat$xintercept_max <- 
    apply(vline_alldat[, grep("id", colnames(vline_alldat), invert = TRUE), drop = FALSE], 1, max) 
  
  if ((!is.null(e_ids)) && (length(e_ids) > 20L)) {
    stop("Unable to generate 'plot_stack' for this many individuals.", call. = FALSE)      
  } else if ((!is.null(e_ids)) && (length(e_ids) > 3L)) {
    warning("'plot_stack' is unlikely to be legible with more than a few individuals.",
            immediate. = TRUE, call. = FALSE)
  }
  
  graph_facet <- if (!is.null(e_ids)) 
    facet_wrap(~ id, scales = "free", nrow = 1) else NULL
  survplot_updated <- survplot + expand_limits(x = 0) + graph_facet + if (vline_found) 
      geom_vline(aes_string(xintercept = "xintercept_max"), vline_alldat, linetype = 2)

  if (!is.null(e_ids)) 
    yplot <- lapply(yplot, function(x) x + facet_wrap(~ id, scales = "free", nrow = 1))
  
  do.call(cowplot::plot_grid, c(yplot, list(survplot_updated), ncol = 1))
}

