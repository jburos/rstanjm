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

# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanjm)

seed <- 123
set.seed(seed)

iter <- 800
chains <- 1
refresh <- iter

if (interactive()) options(mc.cores = parallel::detectCores())

fixef_tol <- 0.05
ranef_tol <- 0.20 


context("stan_jm returns correct estimates for univariate JM, gaussian outcome")
test_that("stan_jm returns correct estimates for univariate JM, gaussian outcome", {
  fit <- stan_jm(formulaLong = y1 ~ t0 + (1 | id),
                 dataLong = simdata_jm_gaus_cvassoc, 
                 formulaEvent = Surv(stop, died) ~ trt,
                 dataEvent = simdata_jm_gaus_cvassoc,
                 time_var = "time", 
                 assoc = "etavalue",
                 iter = iter, refresh = refresh,
                 chains = chains, seed = seed)
  ans.fixef <- list(
    Long1 = c(0, 0.3),
    Event = c(0.1, -0.5, 0.25))

  expect_equal(fixef(fit)$Long1, ans.fixef$Long1, tol = fixef_tol)
  expect_equal(fixef(fit)$Event, ans.fixef$Event, tol = fixef_tol)
  
})
