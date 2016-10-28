# Part of the rstanjm package
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
#
#' PBC datasets used for \pkg{rstanjm} examples
#' 
#' Description of the Mayo Clinic's primary billary chirosis (PBC) 
#' dataset, used for examples in the \pkg{rstanjm} package.
#'
#' @name pbc-datasets
#' @keywords internal
#' @aliases pbcLong pbcSurv pbcLong_subset pbcSurv_subset
#' 
#' @description Longitudinal and time-to-event data for 312 individuals 
#' with primary biliary cirrhosis, who participated in a randomised 
#' placebo controlled trial of D-penicillamine conducted at the Mayo
#' Clinic between 1974 and 1984. The datasets here are derived from the 
#' \code{\link[survival]{pbc}} dataset in the \pkg{survival} 
#' package.    
#'   
#' @format 
#' \describe{
#' The \code{pbcLong} dataset contains both the longitudinal biomarker
#' data and the event time data with multiple rows per individual (one
#' corresponding to each clinic visit). \cr
#' \cr
#' The \code{pbcSurv} dataset contains only one row per individual with 
#' just the event time data. \cr
#' \cr
#' The \code{pbcLong_subset} and \code{pbcSurv_subset} datasets include 
#' a random subset of 40 (of the total 312) individuals.
#' These datasets are small so that they can be used to run examples
#' quickly. \cr
#' \cr
#' The variable definitions are as follows: \cr
#'   \item{\code{age}}{in years}
#'   \item{\code{albumin}}{serum albumin (g/dl)}
#'   \item{\code{alk.phos}}{alkaline phosphotase (U/liter)}
#'   \item{\code{ascites}}{presence of ascites}
#'   \item{\code{ast}}{aspartate aminotransferase, once called SGOT (U/ml)}
#'   \item{\code{bili}}{serum bilirubin (mg/dl)}
#'   \item{\code{logBili}}{logarithm of serum bilirubin}
#'   \item{\code{chol}}{serum cholesterol (mg/dl)}
#'   \item{\code{copper}}{urine copper (ug/day)}
#'   \item{\code{death}}{indicator of death at endpoint} 
#'   \item{\code{edema}}{0 = no edema; 0.5 = untreated or successfully
#'     treated; 1 = edema despite diuretic therapy}
#'   \item{\code{futime}}{time (in days) between registration and the 
#'     earliest of death, transplantion or censoring}
#'   \item{\code{futimeYears}}{time (in years) between registration and  
#'     the earliest of death, transplantion or censoring}
#'   \item{\code{hepato}}{presence of hepatomegaly or enlarged liver}
#'   \item{\code{id}}{numeric ID unique to each individual}
#'   \item{\code{platelet}}{platelet count}
#'   \item{\code{protime}}{standardised blood clotting time}
#'   \item{\code{sex}}{gender (m = male, f = female)}
#'   \item{\code{stage}}{histologic stage of disease (needs biopsy)}
#'   \item{\code{status}}{status at endpoint (0 = censored, 
#'     1 = transplant, 2 = dead)}
#'   \item{\code{trig}}{triglycerides (mg/dl)}
#'   \item{\code{trt}}{binary treatment code (0 = placebo, 1 = 
#'     D-penicillamine)}  
#' }
#' 
NULL
