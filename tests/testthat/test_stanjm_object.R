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

fixef_tol <- 0.10
ranef_tol <- 0.10 


context("Check arguments or slots for stan_jm objects")

# Gaussian models

l1 <- lmer(y1 ~ t0 + (1 | id),
           data = simdata_jm_gaus_cvassoc)

fit1 <- stan_jm(formulaLong = y1 ~ t0 + (1 | id),
                dataLong = simdata_jm_gaus_cvassoc, 
                formulaEvent = Surv(stop, died) ~ trt,
                dataEvent = simdata_jm_gaus_cvassoc.id,
                time_var = "t0", 
                assoc = "etavalue",
                iter = iter, refresh = refresh,
                chains = chains, seed = seed)

fit2 <- stan_jm(formulaLong = y1 ~ t0 + (1 | id),
                dataLong = simdata_jm_gaus_cvassoc, 
                formulaEvent = Surv(start, stop, died) ~ trt,
                dataEvent = simdata_jm_gaus_cvassoc,
                time_var = "t0", 
                assoc = "etavalue",
                iter = iter, refresh = refresh,
                chains = chains, seed = seed)

#--------------------#

test_that("n_subjects argument is right", {
  expect_equal(fit1$n_subjects, ngrps(l1)[["id"]])
  expect_equal(fit2$n_subjects, ngrps(l1)[["id"]])
})

test_that("n_grps argument is right", {
  expect_equal(fit1$n_grps, ngrps(fit1))
  expect_equal(fit2$n_grps, ngrps(fit2))
})







