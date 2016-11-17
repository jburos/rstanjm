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
context("check arguments for stanjm object")


#--------  Models

if (!exists("examplejm")) example(examplejm)

if (!exists("l1")) l1 <- lmer(logBili ~ trt + year + (1 | id), data = pbcLong_subset)

if (!exists("s1")) s1 <- coxph(Surv(futimeYears, death) ~ trt, data = pbcSurv_subset)


#--------  Helper functions

submodel_names <- c("Long1", "Event")
n_submodels <- 2

marker_names <- c("Long1")
n_markers <- 1


#--------  Tests

test_that("coefficients argument returns correct structure", {
  
  expect_equal(length(examplejm$coefficients), n_submodels)
  expect_identical(names(examplejm$coefficients), submodel_names)
  
  # Longitudinal coefficients
  n_lpars <- length(unlist(ranef(l1))) + length(fixef(l1))
  expect_equal(length(examplejm$coefficients$Long1), n_lpars)
  
  # Event coefficients
  n_spars <- length(fixef(s1))
  #expect_equal(length(examplejm$coefficients$Event), n_pars)
  
})


test_that("ses argument returns correct structure", {
  
  expect_equal(length(examplejm$ses), n_submodels)
  expect_identical(names(examplejm$ses), submodel_names)
  
  # Longitudinal coefficients
  n_lpars <- length(unlist(ranef(l1))) + length(fixef(l1))
  expect_equal(length(examplejm$ses$Long1), n_lpars)
  
  # Event coefficients
  n_spars <- length(fixef(s1))
  #expect_equal(length(examplejm$ses$Event), n_pars)
  
})


test_that("fitted.values argument returns correct structure", {

  expect_equal(length(examplejm$fitted.values), n_markers)
  expect_identical(names(examplejm$fitted.values), marker_names)
  expect_equal(length(examplejm$fitted.values$Long1), nobs(l1))
  
})


test_that("linear.predictors argument returns correct structure", {
  
  expect_equal(length(examplejm$linear.predictors), n_markers)
  expect_identical(names(examplejm$linear.predictors), marker_names)
  expect_equal(length(examplejm$linear.predictors$Long1), nobs(l1))
  
})


test_that("residuals argument returns correct structure", {
  
  expect_equal(length(examplejm$residuals), n_markers)
  expect_identical(names(examplejm$residuals), marker_names)
  expect_equal(length(examplejm$residuals$Long1), nobs(l1))
  
})


test_that("n_* arguments return correct values", {
  
  expect_equal(examplejm$n_subjects, ngrps(l1)[["id"]])
  expect_equal(examplejm$n_subjects, ngrps(examplejm)[["id"]])
  expect_equal(examplejm$n_grps, ngrps(l1))
  expect_equal(examplejm$n_grps, ngrps(examplejm))
  expect_equal(examplejm$n_yobs, nobs(l1))
  expect_equal(examplejm$n_events, s1$nevent)

})



