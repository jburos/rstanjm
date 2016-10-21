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
               formulaEvent = Surv(stop, died) ~ trt,
               dataEvent = simdata_jm_gaus_cvassoc,
               time_var = "t0", 
               assoc = "etavalue",
               iter = iter, refresh = refresh,
               chains = chains, seed = seed)


context("stan_jm returns correct attributes and ")
test_that("stan_jm returns correct estimates for univariate JM, gaussian outcome", {

  
  expect_equal(fixef(fit)$Event, ans.fixef$Event, tol = fixef_tol)
  
})

test_that("weights arguments works", {
  idvec0 <- simdata_jm_gaus_cvassoc.id[["id"]]
  idvec1 <- head(idvec0)            # missing IDs
  idvec2 <- rep(idvec0, each = 2)   # repeated IDs
  idvec3 <- c(idvec, 9998, 9999)    # extra IDs not in model
  wts0 <- data.frame(id = idvec0, weights = rep_len(c(1,2), length(idvec0)))
  wts1 <- data.frame(id = idvec1, weights = rep_len(c(1,2), length(idvec1)))
  wts2 <- data.frame(id = idvec2, weights = rep_len(c(1,2), length(idvec2)))
  wts3 <- data.frame(id = idvec0, weights = rep_len(c(1,2), length(idvec0)),
                     junkcol = idvec0)
  wts4 <- data.frame(id = idvec0, weights = rep_len(c("word"), length(idvec0)))
  wts5 <- data.frame(id = idvec0, weights = rep_len(c(NA), length(idvec0)))
  wts6 <- data.frame(id = idvec0, weights = rep_len(c(-1, 1), length(idvec0)))
  wts7 <- data.frame(id = idvec3, weights = rep_len(c(1,2), length(idvec3)))
  
  expect_silent(update(fit1, weights = wts0, iter = 5))
  expect_error(update(fit1, weights = as.matrix(wts0)), "should be a data frame")
  expect_error(update(fit1, weights = wts1), "do not have weights supplied")
  expect_error(update(fit1, weights = wts2), "should only have one row")
  expect_error(update(fit1, weights = wts3), "should be a data frame with two columns")
  expect_error(update(fit1, weights = wts4), "weights supplied must be numeric")
  expect_error(update(fit1, weights = wts5), "weights supplied must be numeric")
  expect_error(update(fit1, weights = wts6), "Negative weights are not allowed")
  expect_silent(update(fit1, weights = wts7, iter = 5)) # ok to supply extra IDs in weights
  
})

