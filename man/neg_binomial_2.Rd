% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neg_binomial_2.R
\name{neg_binomial_2}
\alias{neg_binomial_2}
\title{Negative binomial family function}
\usage{
neg_binomial_2(link = "log")
}
\arguments{
\item{link}{The same as for \code{\link{poisson}}, typically a character
vector of length one among \code{"log"}, \code{"identity"}, and
\code{"sqrt"}.}
}
\value{
An object of class \code{\link[stats]{family}} very similar to
  that of \code{\link[stats]{poisson}} but with a different family name.
}
\description{
Returns the information required to fit a negative binomial generalised
linear mixed model. A call to this function can be passed to the \code{family}
argument of \code{\link{stan_jm}} to specify a negative binomial generalised
linear mixed model for one of the longitudinal submodels. Note that the 
overdispersion parameter \code{theta} is not provided by the user, but
rather it is estimated using the data when the model is fit using 
\code{\link{stan_jm}}.
}
\examples{

}

