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
#' A model for use in the \pkg{rstanjm} examples. 
#' 
#' @name examplejm
#' @format Calling \code{example("examplejm")} will run the model in the 
#'   \strong{Examples} section below, and the resulting stanjm object will 
#'   then be available in the global environment. The example uses
#'   the Mayo Clinic's primary biliary cirrhosis dataset, however, the example
#'   only uses a small subset of the patients so that the model runs
#'   quickly. Type \code{help("pbc-datasets", package = "rstanjm")} for a  
#'   brief description of the data.\cr
#'   \cr
#'   Also, note that by default rstanjm fits models using a single MCMC chain. 
#'   However, the preference is to run multiple MCMC chains to help assess
#'   convergence and also increase the effective sample size. The number of 
#'   chains can be increased by simply specifying the \code{chains} argument
#'   when fitting the model. If you have a multicore CPU with excess RAM, you 
#'   can fit multiple MCMC chains in parallel by setting the \code{cores} 
#'   argument equal to the number of chains being executed. This is preferable 
#'   since then the multiple chains will be run in parallel, thereby saving on 
#'   computation time.
#'   
#' @examples
#' 
#'   examplejm <- 
#'      stan_jm(formulaLong = logBili ~ year + (1 | id), 
#'              dataLong = pbcLong_subset,
#'              formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
#'              dataEvent = pbcSurv_subset,
#'              time_var = "year",
#'              chains = 1, cores = 1, seed = 12345)
#' 
#' 
NULL
