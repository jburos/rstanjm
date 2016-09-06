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
#' 
plot_stack <- function(yplot, survplot) {
    
    y_build <- ggplot_build(yplot)
    y_layout <- y_build$panel$layout
    if (!"id" %in% colnames(y_layout)) ids <- NULL else 
      ids <- y_layout[["id"]]
    
    e_build <- ggplot_build(survplot)
    e_layout <- e_build$panel$layout    
    if (!"id" %in% colnames(e_layout)) e_ids <- NULL else 
      e_ids <- e_layout[["id"]]
    
    if (!identical(ids, e_ids)) 
      stop("The individuals in the 'yplot' and 'survplot' appear to differ. Please ",
           "reestimate the plots using a common 'ids' argument.", call. = FALSE)
    
    L <- length(y_build$data)
    vline_dat <- y_build$data[[L]]
    if (!"xintercept" %in% colnames(vline_dat)) {
      vline <- FALSE
      warning("Could not find vertical line indicating last observation time in the ",
              "plot of the longitudinal trajectory. You may wish to plot the longitudinal ",
              "trajectories again with 'abline = TRUE'.", immediate. = TRUE)     
    } else {
      vline <- TRUE
      vline_dat <- vline_dat[, c("PANEL", "xintercept"), drop = FALSE] 
      if (NROW(y_layout) > 1) {
        panel_id_map <- y_layout[, c("PANEL", "id"), drop = FALSE]
        vline_dat <- merge(vline_dat, panel_id_map, by = "PANEL")        
      }
    }

    if ((!is.null(ids)) && (length(ids) > 20L)) {
      stop("Unable to generate 'plot_stack' for this many individuals.", call. = FALSE)      
    } else if ((!is.null(ids)) && (length(ids) > 3L)) {
      warning("'plot_stack' is unlikely to be legible with more than a few individuals.",
              immediate. = TRUE, call. = FALSE)
    }
    
    graph_facet <- if (!is.null(ids)) 
      facet_wrap(~ id, scales = "free", nrow = 1) else NULL
    survplot_updated <- survplot + expand_limits(x = 0) + graph_facet + if (vline) 
        geom_vline(aes_string(xintercept = "xintercept"), vline_dat, linetype = 2)

    if (!is.null(ids)) 
      y_plot <- yplot + facet_wrap(~ id, scales = "free", nrow = 1)
    
    cowplot::plot_grid(yplot, survplot_updated, ncol = 1)
}

