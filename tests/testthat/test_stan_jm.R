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
# package and then running the code below

library(rstanjm)
context("check stan_jm arguments")


#--------  Models

examplejm1 <- 
  stan_jm(formulaLong = logBili ~ year + (1 | id), 
          dataLong = pbcLong_subset,
          formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
          dataEvent = pbcSurv_subset,
          time_var = "year", iter = 1,
          chains = 1, cores = 1, seed = 12345)

examplejm2 <- 
  stan_jm(formulaLong = logBili ~ year + (year | id), 
          dataLong = pbcLong_subset,
          formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
          dataEvent = pbcSurv_subset,
          time_var = "year", iter = 1,
          chains = 1, cores = 1, seed = 54321)

examplejm3 <- 
  stan_jm(formulaLong = list(
            logBili ~ year + (year | id),
            albumin ~ year + (1 | id)),
          dataLong = pbcLong_subset,
          formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
          dataEvent = pbcSurv_subset,
          time_var = "year", iter = 1,
          chains = 1, cores = 1, seed = 56789)

examplejm4 <- 
  stan_jm(formulaLong = list(
            logBili ~ year + (year | id),
            albumin ~ year + (year | id)),
          dataLong = pbcLong_subset,
          formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
          dataEvent = pbcSurv_subset,
          time_var = "year", iter = 1,
          chains = 1, cores = 1, seed = 98765)


#--------  Tests

test_that("formula argument works", {
  
})


test_that("data argument works", {
  
})


test_that("id_var argument works", {

  # Models with a single grouping factor
    
  expect_output(update(examplejm1, id_var = "id"))
  expect_output(expect_warning(update(examplejm1, id_var = "year"), "are not the same; 'id_var' will be ignored"))
  
  # Models with more than one grouping factor
  
  tmpdat <- pbcLong_subset
  tmpdat$practice <- cut(pbcLong_subset$id, c(0,10,20,30,40))
  expect_error(update(examplejm1, formulaLong. = logBili ~ year + (1 | id) + (1 | practice), dataLong = tmpdat), "'id_var' must be specified")
  expect_error(update(examplejm1, formulaLong. = logBili ~ year + (1 | id) + (1 | practice), dataLong = tmpdat, id_var = "year"), "'id_var' must be included as a grouping factor")
  expect_output(update(examplejm1, formulaLong. = logBili ~ year + (1 | id) + (1 | practice), dataLong = tmpdat, id_var = "id"))
  
})


test_that("family argument works", {
  
})


test_that("assoc argument works", {
  
  # Univariate joint models
  
  expect_output(ret <- update(examplejm2, assoc = NULL))
  expect_output(ret <- update(examplejm2, assoc = "null"))
  expect_output(ret <- update(examplejm2, assoc = "etavalue"))
  expect_output(ret <- update(examplejm2, assoc = "etaslope"))
  expect_output(ret <- update(examplejm2, assoc = "muvalue"))
  expect_output(ret <- update(examplejm2, assoc = "muslope"))
  expect_output(ret <- update(examplejm2, assoc = c("etavalue", "etaslope"))) 
  expect_output(ret <- update(examplejm2, assoc = c("etavalue", "muslope"))) 
  expect_output(ret <- update(examplejm2, assoc = c("muvalue", "etaslope"))) 
  expect_output(ret <- update(examplejm2, assoc = c("muvalue", "muslope")))
  
  expect_error(ret <- update(examplejm2, assoc = c("etavalue", "muvalue")), "cannot be specified together")
  expect_error(ret <- update(examplejm2, assoc = c("etaslope", "muslope")), "cannot be specified together")
  
  expect_output(ret <- update(examplejm2, assoc = "shared_b"))
  expect_output(ret <- update(examplejm2, assoc = "shared_b(1)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_b(2)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_b(1:2)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_b(1,2)"))
  
  expect_output(ret <- update(examplejm2, assoc = "shared_coef"))
  expect_output(ret <- update(examplejm2, assoc = "shared_coef(1)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_coef(2)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_coef(1:2)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_coef(1,2)"))
  
  expect_error(ret <- update(examplejm2, assoc = "shared_b(10)"), "greater than the number of")
  expect_error(ret <- update(examplejm2, assoc = "shared_coef(10)"), "greater than the number of")
  expect_error(ret <- update(examplejm2, assoc = c("shared_b(1)", "shared_coef(1)")), "should not be specified in both")
  expect_error(ret <- update(examplejm2, assoc = c("shared_b", "shared_coef")), "should not be specified in both")
 
  expect_output(ret <- update(examplejm2, assoc = list(NULL)))
  expect_output(ret <- update(examplejm2, assoc = list("null")))
  expect_output(ret <- update(examplejm2, assoc = list("etavalue")))
  expect_output(ret <- update(examplejm2, assoc = list("etaslope")))
  expect_output(ret <- update(examplejm2, assoc = list("muvalue")))
  expect_output(ret <- update(examplejm2, assoc = list("muslope")))
  expect_output(ret <- update(examplejm2, assoc = list(c("etavalue", "etaslope")))) 
  expect_output(ret <- update(examplejm2, assoc = list(c("etavalue", "muslope")))) 
  expect_output(ret <- update(examplejm2, assoc = list(c("muvalue", "etaslope")))) 
  expect_output(ret <- update(examplejm2, assoc = list(c("muvalue", "muslope"))))  
  
  expect_error(ret <- update(examplejm2, assoc = NA), "'assoc' should be") 
  expect_error(ret <- update(examplejm2, assoc = 123), "'assoc' should be") 
  expect_error(ret <- update(examplejm2, assoc = c(1,2,3)), "'assoc' should be") 
  
  expect_error(ret <- update(examplejm2, assoc = c("wrong")), "unsupported association type") 
  expect_error(ret <- update(examplejm2, assoc = list("wrong")), "unsupported association type") 
  
  expect_error(ret <- update(examplejm2, assoc = list(NULL, NULL)), "length equal to the number of") 
  expect_error(ret <- update(examplejm2, assoc = list("etavalue", "etavalue")), "length equal to the number of") 
  expect_error(ret <- update(examplejm2, assoc = list(c("etavalue", "etaslope"), "etavalue")), "length equal to the number of") 

  # Multivariate joint models
  
  expect_output(ret <- update(examplejm3, assoc = "etavalue"))
  expect_output(ret <- update(examplejm3, assoc = "etaslope"))
  expect_output(ret <- update(examplejm3, assoc = "muvalue"))
  expect_output(ret <- update(examplejm3, assoc = "muslope"))
  expect_output(ret <- update(examplejm3, assoc = c("etavalue", "etaslope"))) 
  expect_output(ret <- update(examplejm3, assoc = c("etavalue", "muslope"))) 
  expect_output(ret <- update(examplejm3, assoc = c("muvalue", "etaslope"))) 
  expect_output(ret <- update(examplejm3, assoc = c("muvalue", "muslope")))

  expect_output(ret <- update(examplejm3, assoc = list("etavalue")))
  expect_output(ret <- update(examplejm3, assoc = list("etavalue", "etavalue")))
  expect_output(ret <- update(examplejm3, assoc = list(c("etavalue", "etaslope"), "etavalue")))
  expect_output(ret <- update(examplejm3, assoc = list("etavalue", c("etavalue", "etaslope"))))
  expect_output(ret <- update(examplejm3, assoc = list(c("etavalue", "etaslope"), c("muvalue", "muslope"))))
  
  expect_error(ret <- update(examplejm3, assoc = list("wrong", "etavalue")), "unsupported association type")
  expect_error(ret <- update(examplejm3, assoc = list("null", "etavalue", "etaslope")), "length equal to the number of")
  expect_error(ret <- update(examplejm3, assoc = data.frame("etavalue", "etaslope")), "'assoc' should be") 
  
    
})
  

test_that("base_haz argument works", {

  expect_output(update(examplejm1, base_haz = "weibull"))
  expect_output(update(examplejm1, base_haz = "bs"))
  expect_output(update(examplejm1, base_haz = "piecewise"))

  expect_output(update(examplejm1, base_haz = "bs", df = 5))
  expect_output(update(examplejm1, base_haz = "bs", knots = c(1,3,5)))
  expect_output(update(examplejm1, base_haz = "piecewise", df = 5))
  expect_output(update(examplejm1, base_haz = "piecewise", knots = c(1,3,5)))
  
  expect_output(expect_warning(update(examplejm1, base_haz = "weibull", df = 1), "'df' will be ignored"))
  expect_output(expect_warning(update(examplejm1, base_haz = "weibull", knots = 1), "'knots' will be ignored"))
 
  expect_output(update(examplejm1, base_haz = "piecewise", knots = c(1,3,5)))
  
  
  expect_error(update(examplejm1, base_haz = "bs", df = 1), "'df' must be atleast 3")
  expect_error(update(examplejm1, base_haz = "bs", knots = -1), "'knots' must be non-negative")
  expect_error(update(examplejm1, base_haz = "piecewise", knots = -1), "'knots' must be non-negative")
  expect_error(update(examplejm1, base_haz = "piecewise", knots = c(1,2,50)), "'knots' cannot be greater than the largest event time")
  
})


test_that("quadnodes argument works", {
  
  expect_output(update(examplejm1, quadnodes = 7))
  expect_output(update(examplejm1, quadnodes = 11))
  expect_output(update(examplejm1, quadnodes = 15))
  
  expect_error(update(examplejm1, quadnodes = 1), "'quadnodes' must be either 7, 11 or 15")
  expect_error(update(examplejm1, quadnodes = c(1,2)), "should be a numeric vector of length 1")
  expect_error(update(examplejm1, quadnodes = "wrong"), "should be a numeric vector of length 1")

})


test_that("weights argument works", {
  
  idvec0 <- pbcSurv_subset[["id"]]
  idvec1 <- head(idvec0)            # missing IDs
  idvec2 <- rep(idvec0, each = 2)   # repeated IDs
  idvec3 <- c(idvec0, 9998, 9999)    # extra IDs not in model
  
  wts0 <- data.frame(id = idvec0, weights = rep_len(c(1,2), length(idvec0)))
  wts1 <- data.frame(id = idvec1, weights = rep_len(c(1,2), length(idvec1)))
  wts2 <- data.frame(id = idvec2, weights = rep_len(c(1,2), length(idvec2)))
  wts3 <- data.frame(id = idvec0, weights = rep_len(c(1,2), length(idvec0)),
                     junkcol = idvec0)
  wts4 <- data.frame(id = idvec0, weights = rep_len(c("word"), length(idvec0)))
  wts5 <- data.frame(id = idvec0, weights = rep_len(c(NA), length(idvec0)))
  wts6 <- data.frame(id = idvec0, weights = rep_len(c(-1, 1), length(idvec0)))
  wts7 <- data.frame(id = idvec3, weights = rep_len(c(1,2), length(idvec3)))
  
  expect_output(update(examplejm1, weights = wts0, iter = 5))
  expect_output(update(examplejm1, weights = wts7, iter = 5)) # ok to supply extra IDs in weights
  
  expect_error(update(examplejm1, weights = as.matrix(wts0)), "should be a data frame")
  expect_error(update(examplejm1, weights = wts1), "do not have weights supplied")
  expect_error(update(examplejm1, weights = wts2), "should only have one row")
  expect_error(update(examplejm1, weights = wts3), "should be a data frame with two columns")
  expect_error(update(examplejm1, weights = wts4), "weights supplied must be numeric")
  expect_error(update(examplejm1, weights = wts5), "weights supplied must be numeric")
  expect_error(update(examplejm1, weights = wts6), "Negative weights are not allowed")
  
})


test_that("centre argument works", {
  
  expect_output(update(examplejm1, centreLong = TRUE, centreEvent = TRUE))
  expect_output(update(examplejm1, centreLong = TRUE, centreEvent = FALSE))
  expect_output(update(examplejm1, centreLong = FALSE, centreEvent = TRUE))

})


test_that("init argument works", {
  
  expect_output(update(examplejm1, init = "model_based"))
  expect_output(update(examplejm1, init = "0"))
  expect_output(update(examplejm1, init = 0))
  expect_output(update(examplejm1, init = "random"))
  
})


test_that("prior_PD argument works", {
  
  expect_output(update(examplejm1, prior_PD = TRUE))
  
})

test_that("adapt_delta argument works", {
  
  expect_output(update(examplejm1, adapt_delta = NULL))
  expect_output(update(examplejm1, adapt_delta = 0.8))

})

test_that("max_treedepth argument works", {
  
  expect_output(update(examplejm1, max_treedepth = NULL))
  expect_output(update(examplejm1, max_treedepth = 5))
  expect_output(update(examplejm1, max_treedepth = 5L))
  
})

test_that("QR argument works", {
  
  expect_error(update(examplejm1, QR = TRUE), "not yet implemented")
  
})

