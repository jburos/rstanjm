% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{pairs.stanjm}
\alias{pairs.stanjm}
\title{Pairs method for stanjm objects}
\usage{
\method{pairs}{stanjm}(x, ...)
}
\arguments{
\item{x}{A fitted model object returned by the
\code{\link{stan_jm}} modelling function. See \code{\link{stanjm-object}}.}

\item{...}{Arguments to pass to \code{\link[rstan]{pairs.stanfit}}.}
}
\description{
See \code{\link[rstan]{pairs.stanfit}} for details.
}
\details{
See the Details section in \code{\link[rstan]{pairs.stanfit}}.
}
\examples{
if (!exists("examplejm")) example(examplejm)
pairs(examplejm, pars = c("Long1|(Intercept)", "log-posterior"))

}

