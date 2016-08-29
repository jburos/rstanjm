# Part of the rstanjm package for estimating model parameters
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
#' @name example_jm
#' @format Calling \code{example("example_jm")} will run the model in the 
#'   Examples section, below, and the resulting stanjm object will then be
#'   available in the global environment. The \code{chains} and \code{iter}
#'   arguments are specified to make this example be small in size. In practice,
#'   we recommend that they be left unspecified in order to use the default
#'   values (4 and 2000 respectively) or increased if there are convergence
#'   problems. The \code{cores} argument is optional and on a multicore system,
#'   the user may well want to set that equal to the number of chains being
#'   executed.
#'   
#' @examples
#' example_jm <- 
#'   stan_jm(formulaLong = logBili ~ year + (1 | id), 
#'           dataLong = pbcLong,
#'           formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'           dataEvent = pbcSurv,
#'           time_var = "year",
#'           chains = 1, cores = 1, seed = 12345,
#'           iter = 1000, warmup = 800)
#' example_jm
#' 
NULL
