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

#' Example joint longitudinal and time-to-event model
#' 
#' A model for use in \pkg{rstanjm} examples. 
#' 
#' @name examplejm
#' @format Calling \code{example("examplejm")} will run the model in the 
#'   Examples section, below, and the resulting stanjm object will then be
#'   available in the global environment. The \code{cores} argument is 
#'   optional and on a multicore system, the user may well want to set 
#'   that equal to the number of chains being executed.
#'   
#' @examples
#' examplejm <- 
#'   stan_jm(formulaLong = logBili ~ year + (1 | id), 
#'           dataLong = pbcLong_subset,
#'           formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'           dataEvent = pbcSurv_subset,
#'           time_var = "year",
#'           chains = 1, cores = 1, seed = 12345)
#' 
NULL
