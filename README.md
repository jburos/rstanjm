
<!-- README.md is generated from README.Rmd. Please edit that file -->
rstanjm
=======

[![Travis-CI Build Status](https://travis-ci.org/sambrilleman/rstanjm.svg?branch=master)](https://travis-ci.org/sambrilleman/rstanjm) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/rstanjm)](http://www.r-pkg.org/pkg/rstanjm) [![License](https://img.shields.io/badge/License-GPL%20%28%3E=%203%29-brightgreen.svg)](http://www.gnu.org/licenses/gpl-3.0.html)

**rstanjm** is an R package that allows the user to fit joint (shared parameter) models for longitudinal and time-to-event data under a Bayesian framework. Estimation is carried out using the software [Stan](http://mc-stan.org) (via the **rstan** package).

Both univariate (one longitudinal marker) and multivariate (more than one longitudinal marker) joint models can be estimated. The estimated longitudinal submodel(s) are specified as a (multivariate) generalised linear mixed effects model. Therefore, both continuous and non-continuous (e.g. binary or count data) longitudinal markers can be accomodated through a range of possible link functions and error distributions. Multi-level clustered data (e.g. patients within clinics) can be accomodated in the longitudinal submodel, provided that the individual (e.g. patient) is the lowest level of clustering. If a multivariate joint model is specified, then the dependence between the multiple markers is captured through a shared multivariate normal distribution for the random effects.

The event submodel is specified as a proportional hazards model for which the baseline hazard can be a Weibull distribution, piecewise constant, or approximated using B-splines. Various association structures for linking the longitudinal marker to the risk of the event are possible.

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

The `args = "--preclean"` option is necessary to ensure that the Stan model code is compiled correctly. Note that the package does not currently have any supporting vignettes, hence we can only specify `build_vignettes = FALSE`. Once vignettes become available, then these could also be installed by specifying `build_vignettes = TRUE` (assuming LaTeX is appropriately installed on your PC).

Example
-------

In this section we present two very brief examples showing how the **rstanjm** package's main modelling function can be used to fit either a univariate or multivariate shared parameter joint model. In the examples we use the Mayo Clinic's primary biliary cirrhosis (PBC) data. For a description of the datasets type:

``` r
help("pbc-datasets", package = "rstanjm")
```

First, we fit a simple univariate joint model, with a single normally distributed longitudinal marker, an association structure based on the current value of the linear predictor, and a Weibull baseline hazard (note that we fit this model to a small subset of the PBC data so that the example model runs quickly):

``` r
library(rstanjm)
f1 <- stan_jm(formulaLong = logBili ~ year + (year | id), 
              dataLong = pbcLong_subset,
              formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
              dataEvent = pbcSurv_subset,
              time_var = "year")
```

The fitted model is returned as an object of the S3 class `stanjm`. We have a variety of methods and postestimation functions available for this class, including: `print`, `summary`, `plot`, `fixef`, `ranef`, `coef`, `VarCorr`, `posterior_interval`, `update`, and more. See `?print.stanjm`, `?summary.stanjm`, `?plot` and `?stanjm-methods` for more details.

We can also obtain posterior predictions, either for individuals who were used in the estimation or for new individuals. To obtain the estimated longitudinal trajectory we use the function `posterior_predict`, whilst for the estimated survival function we use `posterior_survfit`. By default these functions will use the observed predictor matrices to generate posterior predictions, however, we can obtain out-of-sample predictions by providing new data via the `newdata`, `newdataLong` and `newdataEvent` arguments. The out-of-sample predictions allow us to marginalise over the distribution of the random effects thereby capturing the estimated between-individual variability. See `?posterior_traj` and `?posterior_survfit` for more details.

Second, we demonstrate how we can fit a more complex shared parameter joint model. In this example we fit a multivariate joint model, with two normally distributed longitudinal markers, an association structure based on the current value and current slope of the linear predictor from the first longitudinal submodel and the random intercept term from the second longitudinal submodel, and a baseline hazard approximated using B-splines. We use a horseshoe shrinkage prior for the three association parameters, and a Student t prior with 5 degrees of freedom for each of the remaining regression coefficients. Note that since this joint model is relatively more complex and contains a larger number of parameters, we need to fit this joint model to the full PBC data, containing 312 individuals, and therefore this model takes a while longer to run:

``` r
mv1 <- stan_jm(
        formulaLong = list(
          logBili ~ year + (year | id), 
          spiders ~ sex + year + (1 | id)),
        dataLong = pbcLong_subset,
        formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
        dataEvent = pbcSurv_subset,
        time_var = "year",
        assoc = list(c("etavalue", "etaslope"), "shared_b(1)"),
        base_haz = "splines",
        priorLong = student_t(df = 5),
        priorEvent = student_t(df = 5),
        priorAssoc = hs())
```

The examples in this README file are intended to be brief. For further details, and a more comprehensive demonstration of the **rstanjm** package see the vignette:

``` r
vignette('rstanjm', package = 'rstanjm')
```

Bug Reports
-----------

If you find any bugs, please report them via email to [Sam Brilleman](mailto:sam.brilleman@monash.edu).

References
----------

1.  Stan Development Team (2015) Stan Modeling Language Users Guide and Reference Manual. <http://mc-stan.org/documentation/>

2.  Taylor JM, Park Y, Ankerst DP, et al. Real-time individual predictions of prostate cancer recurrence using joint models. *Biometrics*. 2013; **69(1)**: 206–213.
