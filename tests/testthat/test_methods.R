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


context("Check methods for stan_jm objects")

# Gaussian models

l1 <- lmer(y1 ~ t0 + (1 | id),
           data = simdata_jm_gaus_cvassoc)

s1 <- coxph(Surv(stop, died) ~ trt, 
            data = simdata_jm_gaus_cvassoc.id)

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

# Poisson models

l3 <- glmer(y1 ~ t0 + (1 | id),
            data = simdata_jm_pois_cvassoc,
            family = "poisson")

fit3 <- update(fit1, family = "poisson",
               dataLong = simdata_jm_pois_cvassoc,
               dataEvent = simdata_jm_pois_cvassoc.id)


    #--------------------#

test_that("confint and posterior_interval", {
  
  expect_error(confint(fit1), regexp = "use posterior_interval")
  
  expect_silent(ci0 <- posterior_interval(fit1, prob = 0.5))
  expect_silent(ci1 <- posterior_interval(fit1, prob = 0.5, 
                                          regex_pars = c("^Long1", "^Event", "^Assoc")))
  expect_silent(ci2 <- posterior_interval(fit1, prob = 0.95, regex_pars = "^Long1"))
  expect_silent(ci3 <- posterior_interval(fit1, prob = 0.95, regex_pars = "Event"))
  expect_silent(ci4 <- posterior_interval(fit1, prob = 0.95, regex_pars = "Assoc"))
  expect_silent(ci5 <- posterior_interval(fit1, prob = 0.95, regex_pars = "b\\["))
  expect_silent(ci6 <- posterior_interval(fit1, prob = 0.8, pars = "Long1|(Intercept)"))
  expect_silent(ci7 <- posterior_interval(fit1, regex_pars = "b\\[",
                                          pars = c("Long1|(Intercept)", "Event|weibull-shape")))
  
  expect_identical(rownames(ci1), rownames(summary(fit1)))
  b_nms <- rstanjm:::b_names(rownames(fit1$stan_summary), value = TRUE)
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
  expect_error(posterior_interval(fit1, prob = c(0.25, 0.75)), regexp = prob_msg)
  expect_error(posterior_interval(fit1, prob = 0), regexp = prob_msg)
  expect_error(posterior_interval(fit1, prob = 1), regexp = prob_msg)
  expect_error(posterior_interval(fit1, prob = 2), regexp = prob_msg)  
})

test_that("ngrps is right", {
  expect_equal(ngrps(fit1), ngrps(l1))
  expect_equal(ngrps(fit2), ngrps(l1))
})

test_that("vcov returns correct structure", {
  # Gaussian model
  col1 <- ncol(vcov(l1)) + ncol(vcov(s1)) + 1  # plus 1 for Weibull intercept
  expect_equal(ncol(vcov(fit1)), col1)
  expect_equal(ncol(vcov(fit2)), col1)
  
  # Poisson model
  col3 <- ncol(vcov(l3)) + ncol(vcov(s1)) + 1 
  expect_equal(ncol(vcov(fit3)), col3)
})

test_that("sigma method works", {
  expect_double <- function(x) expect_type(x, "double")
  
  # need to use :: because sigma is masked by lme4's sigma
  rsigma <- rstanjm::sigma
  expect_double(rsigma(fit1))
  expect_double(rsigma(fit2))
  expect_false(identical(rsigma(fit1), 1))
  expect_false(identical(rsigma(fit2), 1))
  
  # Poisson model -- no sigma
  expect_identical(rsigma(fit3), 1)
})


att_names <- function(object) {
  nms <- names(object)
  att_nms <- names(attributes(object))
  att_nms2 <- lapply(object, function(x) names(attributes(x)))
  c(nms, att_nms, att_nms2)
}
check_att_names <- function(x,y) {
  expect_identical(att_names(x), att_names(y))
}
test_that("VarCorr returns correct structure", {
  vc_fit1 <- VarCorr(fit1); vc_l1 <- VarCorr(l1) 
  vc_fit2 <- VarCorr(fit2)
  vc_fit3 <- VarCorr(fit3); vc_l3 <- VarCorr(l3)
  expect_s3_class(vc_fit1, class(vc_l1))
  expect_s3_class(vc_fit2, class(vc_l1))
  expect_s3_class(vc_fit3, class(vc_l3))
  check_att_names(vc_fit1, vc_l1)
  check_att_names(vc_fit2, vc_l1)
  check_att_names(vc_fit3, vc_l3)
})


check_sizes <- function(x,y) {
  expect_equal(length(x), length(y))
  expect_equal(lapply(x, dim), lapply(y, dim))
}
test_that("ranef returns correct structure", {
  re_fit1 <- ranef(fit1)[[1]]; re_l1 <- ranef(l1)
  re_fit2 <- ranef(fit2)[[1]]
  re_fit3 <- ranef(fit3)[[1]]; re_l1 <- ranef(l3)
  expect_s3_class(re_fit1, class(re_l1))
  expect_s3_class(re_fit2, class(re_l1))
  expect_s3_class(re_fit3, class(re_l3))
  check_att_names(re_fit1, re_l1)
  check_att_names(re_fit2, re_l1)
  check_att_names(re_fit3, re_l3)
  check_sizes(re_fit1, re_l1)
  check_sizes(re_fit2, re_l1)
  check_sizes(re_fit3, re_l3)
})

test_that("fixef returns the right coefs", {
  # Longitudinal coefficients
  expect_identical(names(fixef(fit1)$Long1), names(fixef(l1)))
  expect_identical(names(fixef(fit2)$Long1), names(fixef(l1)))
  expect_identical(names(fixef(fit3)$Long1), names(fixef(l3)))
  
  # Event coefficients
  nms1 <- grep("Intercept|Assoc", names(fixef(fit1)$Event), 
               value = TRUE, invert = TRUE)
  nms2 <- grep("Intercept|Assoc", names(fixef(fit2)$Event), 
               value = TRUE, invert = TRUE)
  nms3 <- grep("Intercept|Assoc", names(fixef(fit3)$Event), 
               value = TRUE, invert = TRUE)
  expect_identical(nms1, names(coef(s1)))
  expect_identical(nms2, names(coef(s1)))
  expect_identical(nms3, names(coef(s1)))
})

test_that("coef returns the right structure", {
  # Longitudinal coefficients
  coef_fit1 <- coef(fit1)$Long1; coef_l1 <- coef(l1)
  coef_fit2 <- coef(fit2)$Long1; coef_l1 <- coef(l1)
  coef_fit3 <- coef(fit3)$Long1; coef_l1 <- coef(l1)
  expect_s3_class(coef_fit1, class(coef_l1))
  expect_s3_class(coef_fit2, class(coef_l1))
  expect_s3_class(coef_fit3, class(coef_l3))
  check_att_names(coef_fit1, coef_l1)
  check_att_names(coef_fit2, coef_l1)
  check_att_names(coef_fit3, coef_l3)
  check_sizes(coef_fit1, coef_l1)
  check_sizes(coef_fit2, coef_l1)
  check_sizes(coef_fit3, coef_l3)
  
  # Event coefficients
  nms1 <- grep("Intercept|Assoc", names(coef(fit1)$Event), 
               value = TRUE, invert = TRUE)
  nms2 <- grep("Intercept|Assoc", names(coef(fit2)$Event), 
               value = TRUE, invert = TRUE)
  nms3 <- grep("Intercept|Assoc", names(coef(fit3)$Event), 
               value = TRUE, invert = TRUE)
  expect_identical(nms1, names(coef(s1)))
  expect_identical(nms2, names(coef(s1)))
  expect_identical(nms3, names(coef(s1)))
})

test_that("coef ok if any 'ranef' missing from 'fixef'", {
  l99 <- update(l1, formula = . ~ (t0 | id))
  fit99 <- update(fit1, formulaLong = . ~ (t0 | id))
  coef_fit99 <- coef(fit99)$Long1; coef_l99 <- coef(l99)
  check_att_names(coef_fit99, coef_l99)
  check_sizes(coef_fit99, coef_l99)
})


