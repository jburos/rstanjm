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

#' Prior distributions and options for fitting Bayesian joint models via Stan
#' 
#' These functions are used to specify priors and related options 
#' when fitting a joint model using the \code{\link{stan_jm}} 
#' function, which is the main modeling function in the \pkg{rstanjm} package.
#' The default priors are intended to be 
#' \emph{weakly informative} in that they provide moderate regularlization and 
#' help stabilize computation. For many applications the defaults will perform 
#' well, but prudent use of more informative priors is encouraged. Uniform prior
#' distributions are possible by setting the prior to \code{NULL} but, 
#' unless the data is very strong,
#' they are not recommended and are \emph{not} non-informative, since they give
#' the same probability mass to implausible values as plausible ones.
#' 
#' @name priors
#' @export 
#'   
#' @param location Prior location. For \code{normal} and \code{student_t} 
#'   (provided that \code{df > 1}) this is the prior mean. For \code{cauchy} 
#'   (which is equivalent to \code{student_t} with \code{df=1}), the mean does 
#'   not exist and \code{location} is the prior median. The default value is 
#'   \eqn{0}.
#' @param scale Prior scale. The default depends on the family (see Details).
#' @param df,df1,df2 Prior degrees of freedom. The default is \eqn{1} for 
#'   \code{student_t}, in which case it is equivalent to \code{cauchy}. For the
#'   hierarchical shrinkage priors (\code{hs} and \code{hs_plus}) the degrees of
#'   freedom parameter(s) default to \eqn{3}.
#'   
#' @details The details depend on the family of the prior being used:
#' \subsection{Student t family}{
#'   Family members:
#'   \itemize{
#'   \item \code{normal(location, scale)}
#'   \item \code{student_t(df, location, scale)}
#'   \item \code{cauchy(location, scale)}
#'   }
#'   
#'   For the prior distribution for the intercept, \code{location}, 
#'   \code{scale}, and \code{df} should be scalars. For the prior for the other
#'   coefficients they can either be vectors of length equal to the number of
#'   coefficients (not including the intercept), or they can be scalars, in 
#'   which case they will be recycled to the appropriate length. As the 
#'   degrees of freedom approaches infinity, the Student t distribution 
#'   approaches the normal distribution and if the degrees of freedom are one,
#'   then the Student t distribution is the Cauchy distribution.
#'   
#'   If \code{scale} is not specified it will default to 10 for the intercept
#'   and 2.5 for the other coefficients, unless the probit link function is
#'   used, in which case these defaults are scaled by a factor of 
#'   \code{dnorm(0)/dlogis(0)}, which is roughly 1.6.
#' }
#' \subsection{Hierarchical shrinkage family}{
#'   Family members:
#'   \itemize{
#'   \item \code{hs(df)}
#'   \item \code{hs_plus(df1, df2)}
#'   }
#'   
#'   The hierarchical shrinkage priors are normal with a mean of zero and a 
#'   standard deviation that is also a random variable. The traditional 
#'   hierarchical shrinkage prior utilizes a standard deviation that is 
#'   distributed half Cauchy with a median of zero and a scale parameter that is
#'   also half Cauchy. This is called the "horseshoe prior". The hierarchical 
#'   shrinkage (\code{hs}) prior in the \pkg{rstanjm} package instead utilizes 
#'   a half Student t distribution for the standard deviation (with 3 degrees of
#'   freedom by default), scaled by a half Cauchy parameter, as described by
#'   Piironen and Vehtari (2015). It is possible to change the \code{df}
#'   argument, the prior degrees of freedom, to obtain less or more shrinkage.
#'   
#'   The hierarhical shrinkpage plus (\code{hs_plus}) prior is a normal with a 
#'   mean of zero and a standard deviation that is distributed as the product of
#'   two independent half Student t parameters (both with 3 degrees of freedom
#'   (\code{df1}, \code{df2}) by default) that are each scaled by the same
#'   square root of a half Cauchy parameter.
#'   
#'   These hierarchical shrinkage priors have very tall modes and very fat 
#'   tails. Consequently, they tend to produce posterior distributions that are
#'   very concentrated near zero, unless the predictor has a strong influence on
#'   the outcome, in which case the prior has little influence. Hierarchical
#'   shrinkage priors often require you to increase the 
#'   \code{\link{adapt_delta}} tuning parameter in order to diminish the number
#'   of divergent transitions. For more details on tuning parameters and
#'   divergent transitions see the Troubleshooting section of the 
#'   \emph{How to Use the rstanarm Package} vignette.
#' }
#' \subsection{Covariance matrices}{
#'   Family members:
#'   \itemize{
#'   \item \code{decov(regularization, concentration, shape, scale)}
#'   }
#'   
#'   Covariance matrices are decomposed into correlation matrices and 
#'   variances. The variances are in turn decomposed into the product of a
#'   simplex vector and the trace of the matrix. Finally, the trace is the
#'   product of the order of the matrix and the square of a scale parameter.
#'   This prior on a covariance matrix is represented by the \code{decov} 
#'   function.
#'   
#'   The prior for a correlation matrix is called LKJ whose density is 
#'   proportional to the determinant of the correlation matrix raised to the 
#'   power of a positive regularization parameter minus one. If
#'   \code{regularization = 1} (the default), then this prior is jointly 
#'   uniform over all correlation matrices of that size. If 
#'   \code{regularization > 1}, then the identity matrix is the mode and in the
#'   unlikely case that \code{regularization < 1}, the identity matrix is the
#'   trough.
#'   
#'   The trace of a covariance matrix is equal to the sum of the variances. We
#'   set the trace equal to the product of the order of the covariance matrix
#'   and the \emph{square} of a positive scale parameter. The particular
#'   variances are set equal to the product of a simplex vector --- which is
#'   non-negative and sums to \eqn{1} --- and the scalar trace. In other words,
#'   each element of the simplex vector represents the proportion of the trace
#'   attributable to the corresponding variable.
#'   
#'   A symmetric Dirichlet prior is used for the simplex vector, which has a 
#'   single (positive) \code{concentration} parameter, which defaults to
#'   \eqn{1} and implies that the prior is jointly uniform over the space of
#'   simplex vectors of that size. If \code{concentration > 1}, then the prior
#'   mode corresponds to all variables having the same (proportion of total)
#'   variance, which can be used to ensure the the posterior variances are not
#'   zero. As the \code{concentration} parameter approaches infinity, this
#'   mode becomes more pronounced. In the unlikely case that 
#'   \code{concentration < 1}, the variances are more polarized.
#'   
#'   If all the variables were multiplied by a number, the trace of their 
#'   covariance matrix would increase by that number squared. Thus, it is 
#'   reasonable to use a scale-invariant prior distribution for the positive
#'   scale parameter, and in this case we utilize a Gamma distribution, whose
#'   \code{shape} and \code{scale} are both \eqn{1} by default, implying a
#'   unit-exponential distribution. Set the \code{shape} hyperparameter to some
#'   value greater than \eqn{1} to ensure that the posterior trace is not zero.
#'   
#'   If \code{regularization}, \code{concentration}, \code{shape} and / or 
#'   \code{scale} are positive scalars, then they are recycled to the 
#'   appropriate length. Otherwise, each can be a positive vector of the 
#'   appropriate length, but the appropriate length depends on the number of 
#'   covariance matrices in the model and their sizes. A one-by-one covariance 
#'   matrix is just a variance and thus does not have \code{regularization} or 
#'   \code{concentration} parameters, but does have \code{shape} and 
#'   \code{scale} parameters for the prior standard deviation of that 
#'   variable.
#' }
#'   
#' @return A named list to be used internally by the model
#'   fitting function \code{\link{stan_jm}}.
#'   
normal <- function(location = 0, scale = NULL) {
  validate_parameter_value(scale)
  nlist(dist = "normal", df = NA, location, scale)
}

#' @rdname priors
#' @export
student_t <- function(df = 1, location = 0, scale = NULL) {
  validate_parameter_value(scale)
  validate_parameter_value(df)
  nlist(dist = "t", df, location, scale)
}

#' @rdname priors
#' @export
cauchy <- function(location = 0, scale = NULL) {
  student_t(df = 1, location = location, scale = scale)
}

#' @rdname priors
#' @export
hs <- function(df = 3) {
  validate_parameter_value(df)
  nlist(dist = "hs", df, location = 0, scale = 1)
}

#' @rdname priors
#' @export
hs_plus <- function(df1 = 3, df2 = 3) {
  validate_parameter_value(df1)
  validate_parameter_value(df2)
  # scale gets used as a second df hyperparameter
  nlist(dist = "hs_plus", df = df1, location = 0, scale = df2)
}

#' @rdname priors
#' @export
#' @param regularization Exponent for an LKJ prior on the correlation matrix in
#'   the \code{decov} prior. The default is \eqn{1}, implying a joint uniform
#'   prior.
#' @param concentration Concentration parameter for a symmetric Dirichlet 
#'   distribution. The defaults is \eqn{1}, implying a joint uniform prior.
#' @param shape Shape parameter for a gamma prior on the scale parameter in the
#'   \code{decov} prior. If \code{shape} and \code{scale} are both \eqn{1} (the
#'   default) then the gamma prior simplifies to the unit-exponential
#'   distribution.
decov <- function(regularization = 1, concentration = 1, 
                  shape = 1, scale = 1) {
  validate_parameter_value(regularization)
  validate_parameter_value(concentration)
  validate_parameter_value(shape)
  validate_parameter_value(scale)
  nlist(dist = "decov", regularization, concentration, shape, scale)
}

#' @rdname priors
#' @export 
#' @param prior_scale_for_dispersion Prior scale for the standard error of the 
#'   regression in Gaussian models, which is given a half-Cauchy prior truncated
#'   at zero.
#' @param min_prior_scale Minimum prior scale for the intercept and 
#'   coefficients.
#' @param scaled A logical scalar, defaulting to \code{TRUE}. If \code{TRUE} the
#'   \code{prior_scale} is further scaled by the range of the predictor if the 
#'   predictor has exactly two unique values and scaled by twice the standard
#'   deviation of the predictor if it has more than two unique values.
#'
priorLong_options <- function(prior_scale_for_dispersion = 5, 
                              min_prior_scale = 1e-12, 
                              scaled = TRUE) {
  validate_parameter_value(prior_scale_for_dispersion)
  validate_parameter_value(min_prior_scale)
  nlist(scaled, min_prior_scale, prior_scale_for_dispersion)
}

#' @rdname priors
#' @export 
#' @param prior_scale_for_basehaz Usage depends on the baseline hazard 
#'   specified in the \code{base_haz} argument of the  
#'   \code{\link{stan_jm}} call. 
#'   If \code{base_haz = "weibull"} then this argument specifies the
#'   prior scale for the shape parameter of the Weibull distribution, 
#'   which is given a half-Cauchy truncated at zero.
#'   If \code{base_haz = "splines"} then this argument specifies the
#'   prior scale(s) for the spline regression coefficients, 
#'   which are given normal distributions with location 0.
#'   If \code{base_haz = "piecewise"} then this argument specifies the
#'   prior scale for the parameters corresponding to the log baseline 
#'   hazard within each interval, which are given normal distributions
#'   with location 0.
#'  
priorEvent_options <- function(prior_scale_for_basehaz = 5,
                               min_prior_scale = 1e-12, 
                               scaled = TRUE) {
  validate_parameter_value(prior_scale_for_basehaz)
  validate_parameter_value(min_prior_scale)
  nlist(scaled, min_prior_scale, prior_scale_for_basehaz)
}

#' @rdname priors
#' @export 
#'
priorAssoc_options <- function(min_prior_scale = 1e-12) {
  validate_parameter_value(min_prior_scale)
  nlist(min_prior_scale)
}


