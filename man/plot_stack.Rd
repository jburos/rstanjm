% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_stack.R
\name{plot_stack}
\alias{plot_stack}
\title{Combine the plots of the estimated longitudinal trajectory and survival function}
\usage{
plot_stack(yplot, survplot)
}
\arguments{
\item{yplot}{An object of class \code{plot.predict.stanjm}, returned by a 
call to the generic plotting method for objects of class 
\code{predict.stanjm}. If there is more than one longitudinal outcome,
then a list of such objects can be provided.}

\item{survplot}{An object of class \code{plot.survfit.stanjm}, returned by a 
call to the generic plotting method for objects of class 
\code{survfit.stanjm}.}
}
\value{
A single \code{ggplot} object that includes plots of the estimated 
  subject-specific longitudinal trajectories stacked on top of the 
  associated subject-specific survival curve.
}
\description{
This function takes arguments containing the plots of the estimated  
subject-specific longitudinal trajectory (or trajectories if a multivariate  
joint model was estimated) and the plot of the estimated subject-specific 
survival function and combines them into a single figure. This is most
easily understood by running the \strong{Examples} below.
}
\examples{
\donttest{
  if (!exists("examplejm")) example(examplejm)
  ps1 <- posterior_survfit(examplejm, ids = c(7,13,16))
  pt1 <- posterior_predict(examplejm, ids = c(7,13,16),
                           interpolate = TRUE, extrapolate = TRUE)
  plot_surv <- plot(ps1) 
  plot_traj <- plot(pt1, vline = TRUE, plot_observed = TRUE)
  plot_stack(plot_traj, plot_surv)
}
 
}
\seealso{
\code{\link{plot.predict.stanjm}}, \code{\link{plot.survfit.stanjm}},
  \code{\link{posterior_predict}}, \code{\link{posterior_survfit}}
}

