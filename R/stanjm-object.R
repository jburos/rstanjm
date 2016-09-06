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

#' Fitted model object in the \pkg{rstanjm} package
#'
#' @name stanjm-object 
#'
#' @description The main model fitting function in the \pkg{rstanjm} package 
#'   (see \code{\link{stan_jm}}) returns an object of class \code{stanjm}, which 
#'   is a list containing the components described below.  
#'   
#' @return The following components must be included in a \code{stanjm} object.
#'   In most cases the components return a named list, with each element of the list 
#'   containing the described item (for example, coefficients, ses, residuals, etc)
#'   for one of the longitudinal submodels or the event submodel. \cr
#'
#' \describe{
#'   \item{\code{coefficients}}{
#'   Point estimates for the longitudinal and event submodels.
#'   As described in \code{\link{print.stanjm}}.
#'   }
#'   \item{\code{ses}}{
#'   Standard errors for the longitudinal and event submodels.
#'   Based on \code{\link[stats]{mad}}, as described in
#'   \code{\link{print.stanjm}}.
#'   }
#'   \item{\code{residuals}}{
#'   Residuals for the longitudinal submodel(s) only. Of type \code{'response'}.
#'   }
#'   \item{\code{fitted.values}}{
#'   Fitted mean values for the longitudinal submodel(s) only, based on the
#'   estimated parameters in \code{coefficients}. The fitted mean values are based
#'   on both the fixed and random effects. For models using a non-identity link
#'   function the linear predictors are transformed by the inverse link.
#'   }
#'   \item{\code{linear.predictors}}{
#'   Linear fit on the link scale, for the longitudinal submodel(s) only. 
#'   For linear models this is the same as \code{fitted.values}.
#'   }
#'   \item{\code{covmat}}{
#'   Variance-covariance matrix for the coefficients (both fixed and random),
#'   evaluated separately for each of the longitudinal and event submodels.
#'   Calculation of the variance-covariance matrix is based on draws from the
#'   posterior distribution.
#'   }
#'   \item{\code{n_markers}}{
#'   The number of longitudinal markers (submodels).
#'   }
#'   \item{\code{n_subjects}}{
#'   The number of individuals in the fitted model.
#'   }
#'   \item{\code{n_grps}}{
#'   The number of levels for each grouping factor (will be equal to 
#'   \code{n_subjects} if the individual-level is the only clustering level).
#'   }
#'   \item{\code{n_yobs}}{
#'   The number of observations for each longitudinal submodel.
#'   }
#'   \item{\code{n_events}}{
#'   The number of non-censored event times.
#'   }   
#'   \item{\code{id_var,time_var}}{
#'   The names of the variables distinguishing between individuals, and 
#'   representing time, in the longitudinal submodel.
#'   } 
#'   \item{\code{cnms}}{
#'   A list of column names of the random effects according to the grouping 
#'   factors, collapsed across all longitudinal submodels. 
#'   See \code{\link[lme4]{mkReTrms}}. These can be obtained separately
#'   for each longitudinal submodel through the components 
#'   \code{object$glmod[[1]]@cnms}, \code{object$glmod[[2]]@cnms}, etc.
#'   }
#'   \item{\code{x}}{
#'   The design matrix (both fixed and random parts) for each of the 
#'   longitudinal submodels, evaluated at the observation times.
#'   }
#'   \item{\code{xq}}{
#'   The design matrices for each of the longitudinal and event submodels 
#'   evaluated at the event/censoring times and the quadrature points. 
#'   The first \code{n_subjects} rows correpond to the event/censoring times,
#'   the second \code{n_subjects} rows correspond to the first quadrature point,
#'   the third \code{n_subjects} rows correspond to the second quadrature point, 
#'   and so on. The design matrices for the longitudinal submodels
#'   include both the fixed and random parts.
#'   }
#'   \item{\code{y}}{
#'   The response for the longitudinal submodel(s).
#'   }  
#'   \item{\code{fr}}{
#'   The model frame (see \code{\link[lme4]{model.frame.merMod}} or 
#'   \code{\link[survival]{model.frame.coxph}}) for each of the longitudinal and event 
#'   submodels, evaluated at the observation times in the original data  
#'   (not evaluated at the quadrature points necessary for fitting the joint 
#'   model). This component also includes the addition of the \code{id_var} and 
#'   \code{time_var} variables in the model frame, even if they weren't 
#'   present in the original model formula.
#'   } 
#'   \item{\code{eventtime,status}}{
#'   The event (or censoring) times and the status indicator for each individual.
#'   }   
#'   \item{\code{dataLong,dataEvent}}{
#'   The \code{dataLong} and \code{dataEvent} arguments.
#'   }
#'   \item{\code{call}}{
#'   The matched call.
#'   }
#'   \item{\code{formula}}{
#'   The model \code{\link[stats]{formula}} for each of the longitudinal
#'   and event submodels.
#'   }
#'   \item{\code{family}}{
#'   The \code{\link[stats]{family}} object used for the longitudinal
#'   submodel(s).
#'   }
#'   \item{\code{base_haz}}{
#'   A list containing information about the baseline hazard.
#'   }
#'   \item{\code{assoc}}{
#'   A list containing information about the association structure for the
#'   joint model.
#'   }
#'   \item{\code{quadnodes}}{
#'   The number of Gauss-Kronrod quadrature nodes used to evaluate the 
#'   cumulative hazard in the joint likelihood function.
#'   }
#'   \item{\code{quadpoints}}{
#'   A list containing the quadrature points for each individual, based on
#'   the event (or censoring) times and the number of \code{quadnodes}.
#'   }   
#'   \item{\code{prior.info}}{
#'   A list with information about the prior distributions used.
#'   }    
#'   \item{\code{algorithm}}{
#'   The estimation method used. Will be \code{"sampling"} for models
#'   estimated using \code{stan_jm}.
#'   }
#'   \item{\code{stanfit,stan_summary}}{
#'   The object of \code{\link[rstan]{stanfit-class}} returned by RStan and a
#'   matrix of various summary statistics from the stanfit object.
#'   }
#'   \item{\code{glmod}}{
#'   The fitted model object(s) returned by a call(s) to \code{\link[lme4]{lmer}}
#'   or \code{\link[lme4]{glmer}} when estimating the separate longitudinal model(s)
#'   prior to fitting the joint model.
#'   }
#'   \item{\code{coxmod}}{
#'   The fitted model object returned by a call to \code{\link[survival]{coxph}}
#'   when estimating the separate time-to-event model prior to fitting the joint model.
#'   This time-to-event model is evaluated using the observed data in \code{dataEvent};
#'   it is not evaluated at the quadrature points necessary for fitting the joint model.
#'   }   
#'}
#'
#' @seealso Objects of this class have a number of generic 
#'   methods described in \code{\link{stanjm-methods}},
#'   as well as \code{\link{print.stanjm}}, \code{\link{summary.stanjm}},
#'   \code{\link{as.matrix.stanjm}}, as well as the non-generic functions
#'   \code{\link{posterior_traj}}, \code{\link{posterior_survfit}}, 
#'   \code{\link{posterior_predict}}, \code{\link{posterior_interval}},
#'   \code{\link{pp_check}} and \code{\link{ps_check}}.
#'   
NULL