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
context("check stanjm methods")


#--------  Models

if (!exists("examplejm")) example(examplejm)

if (!exists("l1")) l1 <- lmer(logBili ~ trt + year + (1 | id), data = pbcLong_subset)

if (!exists("s1")) s1 <- coxph(Surv(futimeYears, death) ~ trt, data = pbcSurv_subset)

if (!exists("poissonjm")) poissonjm <- stan_jm(
  formulaLong = y1pois ~ t0 + (1 | id),
  dataLong = simdata_jm_cvassoc,
  formulaEvent = Surv(t0, t, d) ~ trt,
  dataEvent = simdata_jm_cvassoc,
  time_var = "t0", iter = 2)


#--------  Helper functions

att_names <- function(object) {
  nms <- names(object)
  att_nms <- names(attributes(object))
  att_nms2 <- lapply(object, function(x) names(attributes(x)))
  c(nms, att_nms, att_nms2)
}

check_att_names <- function(x,y) {
  expect_identical(att_names(x), att_names(y))
}

check_sizes <- function(x,y) {
  expect_equal(length(x), length(y))
  expect_equal(lapply(x, dim), lapply(y, dim))
}


#--------  Tests

test_that("coef returns correct structure", {
  
  # Longitudinal coefficients
  coef_fit1 <- coef(examplejm)$Long1
  coef_l1 <- coef(l1)
  expect_s3_class(coef_fit1, class(coef_l1))
  check_att_names(coef_fit1, coef_l1)
  check_sizes(coef_fit1, coef_l1)
  
  # Event coefficients
  nms1 <- grep("Intercept|Assoc", names(coef(examplejm)$Event), value = TRUE, invert = TRUE)
  expect_identical(nms1, names(coef(s1)))
  
})


test_that("coef returns correct structure even if no corresponding fixed effect component for a random effect", {
  
  l99 <- update(l1, formula = . ~ (t0 | id))
  fit99 <- update(examplejm, formulaLong = . ~ (t0 | id), iter = 2)
  coef_fit99 <- coef(fit99)$Long1
  coef_l99 <- coef(l99)
  check_att_names(coef_fit99, coef_l99)
  check_sizes(coef_fit99, coef_l99)
  
})


test_that("fitted returns correct structure", {
  
  expect_equal(length(fitted), examplejm$n_markers)
  expect_equal(length(fitted(examplejm)$Long1), nobs(l1))
  
})


test_that("confint and posterior_interval", {
  
  expect_error(confint(examplejm), regexp = "use posterior_interval")
  
  expect_silent(ci0 <- posterior_interval(examplejm, prob = 0.5))
  expect_silent(ci1 <- posterior_interval(examplejm, prob = 0.5, 
                                          regex_pars = c("^Long1", "^Event", "^Assoc")))
  expect_silent(ci2 <- posterior_interval(examplejm, prob = 0.95, regex_pars = "^Long1"))
  expect_silent(ci3 <- posterior_interval(examplejm, prob = 0.95, regex_pars = "Event"))
  expect_silent(ci4 <- posterior_interval(examplejm, prob = 0.95, regex_pars = "Assoc"))
  expect_silent(ci5 <- posterior_interval(examplejm, prob = 0.95, regex_pars = "b\\["))
  expect_silent(ci6 <- posterior_interval(examplejm, prob = 0.8, pars = "Long1|(Intercept)"))
  expect_silent(ci7 <- posterior_interval(examplejm, regex_pars = "b\\[",
                                          pars = c("Long1|(Intercept)", "Event|weibull-shape")))
  
  expect_identical(rownames(ci1), rownames(summary(examplejm)))
  b_nms <- rstanjm:::b_names(rownames(examplejm$stan_summary), value = TRUE)
  expect_identical(rownames(ci5), b_nms[-length(b_nms)])
  expect_identical(rownames(ci6), c("Long1|(Intercept)"))
  expect_identical(rownames(ci7), c("Long1|(Intercept)", "Event|weibull-shape",
                                    b_nms[-length(b_nms)]))

  expect_identical(colnames(ci0), c("25%", "75%"))
  expect_identical(colnames(ci2), c("2.5%", "97.5%"))
  expect_identical(colnames(ci6), c("10%", "90%"))
  expect_identical(colnames(ci7), c("5%", "95%"))
  
  expect_error(posterior_interval(lm(mpg ~ wt, data = mtcars)), 
               regexp = "not a stanjm object")
  
  prob_msg <- "'prob' should be a single number greater than 0 and less than 1."
  expect_error(posterior_interval(examplejm, prob = c(0.25, 0.75)), regexp = prob_msg)
  expect_error(posterior_interval(examplejm, prob = 0), regexp = prob_msg)
  expect_error(posterior_interval(examplejm, prob = 1), regexp = prob_msg)
  expect_error(posterior_interval(examplejm, prob = 2), regexp = prob_msg) 
  
})


test_that("ngrps returns correct number", {
  
  expect_equal(ngrps(examplejm), ngrps(l1))
  expect_equal(ngrps(examplejm), ngrps(l1))

})


test_that("vcov returns correct structure", {

  col1 <- ncol(vcov(l1)) + ncol(vcov(s1)) + 1  # plus 1 for Weibull intercept
  expect_equal(ncol(vcov(fit1)), col1)
  expect_equal(ncol(vcov(fit2)), col1)

})


test_that("sigma method works", {
  
  # need to use :: because sigma may be masked by lme4's sigma
  rsigma <- rstanjm::sigma
  
  # Gaussian model -- estimated sigma
  expect_type(rsigma(examplejm), "double")
  expect_false(identical(rsigma(examplejm), 1))

  # Poisson model -- no sigma
  expect_identical(rsigma(poissonjm), 1)
  
})


test_that("VarCorr returns correct structure", {
  
  vc_jm <- VarCorr(examplejm) 
  vc_l1 <- VarCorr(l1) 
  expect_s3_class(vc_jm, class(vc_l1))
  check_att_names(vc_jm, vc_l1)
  
})


test_that("ranef returns correct structure", {
  
  re_jm <- ranef(examplejm)$Long1
  re_l1 <- ranef(l1)
  expect_s3_class(re_jm, class(re_l1))
  check_att_names(re_jm, re_l1)
  check_sizes(re_jm, re_l1)

})


test_that("fixef returns correct structure", {
  
  # Longitudinal coefficients
  expect_identical(names(fixef(examplejm)$Long1), names(fixef(l1)))

  # Event coefficients
  nms1 <- grep("Intercept|Assoc", names(fixef(examplejm)$Event), value = TRUE, invert = TRUE)
  expect_identical(nms1, names(coef(s1)))

})



