
<!-- README.md is generated from README.Rmd. Please edit that file -->
rstanjm
=======

[![Travis-CI Build Status](https://travis-ci.org/sambrilleman/rstanjm.svg?branch=master)](https://travis-ci.org/sambrilleman/rstanjm) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/rstanjm)](http://www.r-pkg.org/pkg/rstanjm) [![License](https://img.shields.io/badge/License-GPL%20%28%3E=%203%29-brightgreen.svg)](http://www.gnu.org/licenses/gpl-3.0.html)

**rstanjm** is an R package which allows the user to fit joint (shared parameter) models for longitudinal and time-to-event data under a Bayesian framework. Estimation is carried out using the software [Stan](http://mc-stan.org) (via the **rstan** package).

Both univariate (one longitudinal marker) and multivariate (more than one longitudinal marker) joint models can be estimated. Both continuous and non-continuous (e.g. binary or count data) longitudinal markers can be accomodated through a range of possible link functions and error distributions. Multi-level clustered data (e.g. patients within clinics) can be accomodated in the longitudinal submodel, provided that the individual (e.g. patient) is the lowest level of clustering.

A range of prior distributions are available for the regression coefficients (based on those available in the **rstanarm** package). Post-estimation functions, for example posterior predictions for the longitudinal and survival outcomes, are also provided.

**Note:** Please note that the version available on GitHub is the most up-to-date *development* version of the package. A stable version of the package will be available from CRAN once it is released.

Getting Started
---------------

### Prerequisites

The **rstanjm** package requires the C++ toolchain, **rstan** package, and **rstanarm** package to all be appropriately installed. A convenient shortcut to try would be to install the **rstanarm** package first, directly from CRAN (since **rstanarm** also requires the C++ toolchain and includes **rstan** as a dependency). This can be attempted by executing the following code from within your R session:

``` r
install.packages("rstanarm")
```

If that fails, then you will need to install the C++ toolchain and **rstan** package by following these [instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

### Installation

Once **rstan** is successfully installed, you can install **rstanjm** directly from GitHub using the **devtools** package. To do this you should first check you have devtools installed by executing the following commands from within your R session:

``` r
if (!require(devtools)) {
  install.packages("devtools")
}
```

Then execute the following commands to install **rstanjm**:

``` r
library(devtools)
install_github("sambrilleman/rstanjm", args = "--preclean", local = FALSE, build_vignettes = FALSE)
```

The `args = "--preclean"` option is necessary to ensure that the Stan model code is compiled correctly. Note that the package does not currently have any supporting vignettes, hence we can only specify `build_vignettes = FALSE`. Once vignettes become available, then these could also be installed by specifying `build_vignettes = TRUE` (assuming LaTeX is appropriately installed).

Example
=======

Model fitting
-------------

In this section we present some examples based on a small subset of the Mayo Clinic's primary biliary cirrhosis (PBC) data. For a description of the dataset type:

``` r
help("pbc-datasets", package = "rstanjm")
```

First, we fit a simple univariate joint model, with a single normally distributed longitudinal marker, an association structure based on the current value of the linear predictor, and Weibull baseline hazard.

``` r
library(rstanjm)
f1 <- stan_jm(formulaLong = logBili ~ year + (year | id), 
              dataLong = pbcLong_subset,
              formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
              dataEvent = pbcSurv_subset,
              time_var = "year")
```

Next, we fit a multivariate joint model, with two normally distributed longitudinal markers, an association structure based on the current value and current slope of the linear predictor from the first longitudinal submodel and the random intercept (including the fixed component) from the second longitudinal submodel, and a baseline hazard approximated using B-splines. We use a horseshoe shrinkage prior for the three association parameters, and Student t distribution priors with 5 degrees of freedom for the remaining regression coefficients. Note that since this joint model is relatively more complex and contains a larger number of parameters, we need to fit this joint model to the full PBC data containing 312 individuals and therefore this model takes a little longer to run (~ 10 to 20min).

``` r
mv1 <- stan_jm(
        formulaLong = list(
          logBili ~ year + (year | id), 
          albumin ~ sex + year + (1 | id)),
        dataLong = pbcLong_subset,
        formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
        dataEvent = pbcSurv_subset,
        time_var = "year",
        assoc = list(c("etavalue", "etaslope"), "shared_coef(1)"),
        base_haz = "splines",
        priorLong = student_t(df = 5),
        priorEvent = student_t(df = 5),
        priorAssoc = hs())
```

The fitted model is returned as an object of the S3 class `stanjm`. We have a variety of methods and postestimation functions available for this class, including: `print`, `summary`, `plot`, `fixef`, `ranef`, `coef`, `VarCorr`, `posterior_interval`, `update`, and more. See `?print.stanjm`, `?summary.stanjm`, `?plot` and `?stanjm-methods` to get more details.

In-sample predictions
---------------------

We can obtain posterior predictions (either for individuals who were used in the estimation, or for new individuals). For the longitudinal submodel we can use the function `posterior_traj`, which by default will use the observed predictor matrix to generate posterior predictions for the longitudinal submodel:

``` r
pt <- posterior_traj(f1)
head(pt)
```

To easily interpolate between observation times, so that we obtain predictions at equally spaced times across the entire observation period for each individual we can simply set the argument `interpolate = TRUE`. Similarly if we want to obtain predictions at time points beyond the last observation time we can simply set the argument `extrapolate = TRUE`. We can then plot those trajectories for a subset of individuals. The code required is as follows:

``` r
pt <- posterior_traj(f1, interpolate = TRUE, extrapolate = TRUE)
pt_plot <- plot(pt, ids = 1:3)
pt_plot
```

Similarly, we can obtain posterior predictions for the survival function. Here we will generate the estimated survival function for the same three individuals, and then plot estimated survival function alongside the longitudinal predictions using the function `plot_stack`:

``` r
ps <- posterior_survfit(f1)
ps_plot <- plot(ps, ids = 1:3)
plot_stack(pt_plot, ps_plot)
```

Out-of-sample predictions
-------------------------

If we have new data for an individual who was not included in the estimation sample, then we can easily generate posterior predictions for that individual by including the `newdata` argument as follows:

``` r
df <- pbcLong[pbcLong$id == 100, ]
pt_new <- posterior_traj(f1, newdata = df, int = TRUE, ext = TRUE)
plot(pt_new)
```

And the estimated survival function can be obtained using:

``` r
df_surv <- pbcSurv[pbcSurv$id == 100, ]
ps_new <- posterior_survfit(f1, newdataLong = df, newdataEvent = df_surv)
plot(ps_new)
```

Note however that the out-of-sample predictions that are obtained from using these functions are based on draws of the random effects distribution conditional on the data used to fit the model, but *not* conditional on the new data. That is, they are fitted values evaluated at the set of new covariate values provided, but *not* dynamic predictions (see for example Taylor et al. (2013)).

Model checking
--------------

We can also check the fit of the longitudinal and event submodels by comparing posterior predictions to the observed data. This can be easily done using the `pp_check` and `ps_check` functions. For the unvariate joint model:

``` r
pp_check(f1)
ps_check(f1)
```

or for the multivariate joint model:

``` r
pp_check(mv1, m = 1)  # check fit for first longitudinal submodel
pp_check(mv1, m = 2)  # check fit for second longitudinal submodel
ps_check(mv1)
```

We can also use the diagnostic and inference plots available through the generic `plot` method for `stanjm` objects (see `?plot.stanjm` for details). For example:

``` r
plot(f1, plotfun = "trace")
```

Bug Reports
===========

If you find any bugs, please report them via email to [Sam Brilleman](mailto:sam.brilleman@monash.edu).

References
==========

1.  Stan Development Team (2015) Stan Modeling Language Users Guide and Reference Manual. <http://mc-stan.org/documentation/>

2.  Taylor JM, Park Y, Ankerst DP, et al. Real-time individual predictions of prostate cancer recurrence using joint models. *Biometrics*. 2013; **69(1)**: 206–213.
