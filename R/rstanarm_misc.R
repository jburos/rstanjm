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


# The following functions are taken verbatim from the rstanarm
# package. They are included here, so that they can be used
# internally within the rstanjm package, without needing to
# extract them from the rstanarm namespace using ':::'

#----------  Functions from: misc.R ----------# 

# If a is NULL (and Inf, respectively) return b, otherwise just return a
# @param a,b Objects
`%ORifNULL%` <- function(a, b) {
  if (is.null(a)) b else a
}
`%ORifINF%` <- function(a, b) {
  if (a == Inf) b else a
}

# Combine pars and regex_pars
#
# @param x stanreg object
# @param pars Character vector of parameter names
# @param regex_pars Character vector of patterns
collect_pars <- function(x, pars = NULL, regex_pars = NULL) {
  if (is.null(pars) && is.null(regex_pars)) 
    return(NULL)
  if (!is.null(pars)) 
    pars[pars == "varying"] <- "b"
  if (!is.null(regex_pars)) 
    pars <- c(pars, grep_for_pars(x, regex_pars))
  unique(pars)
}

# Test if an object is a stanreg object
#
# @param x The object to test. 
is.stanreg <- function(x) inherits(x, "stanreg")

# Throw error if object isn't a stanreg object
# 
# @param x The object to test.
validate_stanreg_object <- function(x, call. = FALSE) {
  if (!is.stanreg(x))
    stop("Object is not a stanreg object.", call. = call.) 
}

# Get the posterior sample size
#
# @param x A stanreg object
# @return NULL if used.optimizing(x), otherwise the posterior sample size
posterior_sample_size <- function(x) {
  validate_stanreg_object(x)
  if (used.optimizing(x)) 
    return(NULL)
  pss <- x$stanfit@sim$n_save
  if (used.variational(x))
    return(pss)
  sum(pss - x$stanfit@sim$warmup2)
}

# Test for a given estimation method
#
# @param x A stanreg object.
used.optimizing <- function(x) {
  x$algorithm == "optimizing"
}
used.sampling <- function(x) {
  x$algorithm == "sampling"
}
used.variational <- function(x) {
  x$algorithm %in% c("meanfield", "fullrank")
}

# Test if stanreg object used stan_(g)lmer
#
# @param x A stanreg object.
is.mer <- function(x) {
  check1 <- inherits(x, "lmerMod")
  check2 <- !is.null(x$glmod)
  if (check1 && !check2) {
    stop("Bug found. 'x' has class 'lmerMod' but no 'glmod' component.")
  } else if (!check1 && check2) {
    stop("Bug found. 'x' has 'glmod' component but not class 'lmerMod'.")
  }
  isTRUE(check1 && check2)
}

# Consistent error message to use when something is only available for 
# models fit using MCMC
#
# @param what An optional message to prepend to the default message.
STOP_sampling_only <- function(what) {
  msg <- "only available for models fit using MCMC (algorithm='sampling')."
  if (!missing(what)) 
    msg <- paste(what, msg)
  stop(msg, call. = FALSE)
}

# Consistent error message to use when something is only available for models
# fit using MCMC or VB but not optimization
# 
# @param what An optional message to prepend to the default message.
STOP_not_optimizing <- function(what) {
  msg <- "not available for models fit using algorithm='optimizing'."
  if (!missing(what)) 
    msg <- paste(what, msg)
  stop(msg, call. = FALSE)
}


#----------  Functions from: print-and-summary.R ----------# 

.printfr <- function(x, digits, ...) {
  print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
}
.median_and_madsd <- function(x) {
  cbind(Median = apply(x, 2, median), MAD_SD = apply(x, 2, mad))
}

# Regex parameter selection
#
# @param x stanreg object
# @param regex_pars Character vector of patterns
grep_for_pars <- function(x, regex_pars) {
  validate_stanreg_object(x)
  if (used.optimizing(x)) {
    warning("'regex_pars' ignored for models fit using algorithm='optimizing'.",
            call. = FALSE)
    return(NULL)
  }
  stopifnot(is.character(regex_pars))
  out <- unlist(lapply(seq_along(regex_pars), function(j) {
    grep(regex_pars[j], rownames(x$stan_summary), value = TRUE) 
  }))
  if (!length(out))
    stop("No matches for 'regex_pars'.", call. = FALSE)
  
  return(out)
}


#----------  Functions from: plots.R ----------# 

# Check for valid parameters
# @param x stanreg object
# @param pars user specified character vector
check_plotting_pars <- function(x, pars) {
  if (used.optimizing(x)) {
    allpars <- c("alpha", "beta", rownames(x$stan_summary))
  } else {
    sim <- x$stanfit@sim
    allpars <- c(sim$pars_oi, sim$fnames_oi)
  }
  m <- which(match(pars, allpars, nomatch = 0) == 0)
  if (length(m) > 0) 
    stop("No parameter ", paste(pars[m], collapse = ', '), 
         call. = FALSE) 
  return(unique(pars))
}

# Select the correct plotting function
# @param x stanreg object
# @param plotfun user specified plotfun argument (can be missing)
set_plotting_fun <- function(x, plotfun = NULL) {
  .plotters <- function(x) paste0("stan_", x)
  
  if (used.optimizing(x)) {
    if (!is.null(plotfun)) {
      stop("'plotfun' should not be specified for models fit using ",
           "algorithm='optimizing'.", call. = FALSE)
    } else {
      return("stan_plot_opt")
    }
  } else if (is.null(plotfun)) {
    plotfun <- "stan_plot"
  }
  
  samp_only <- c("ac", "diag", "rhat", "ess", "mcse", "par")
  plotters <- .plotters(c("plot", "trace", "scat", "hist", "dens", samp_only))
  funname <- grep(paste0(plotfun, "$"), plotters, value = TRUE)
  if (used.variational(x) && funname %in% .plotters(samp_only))
    STOP_sampling_only(funname)
  fun <- try(getExportedValue("rstan", funname), silent = TRUE)
  if (inherits(fun, "try-error")) 
    stop("Plotting function not found. See ?rstanarm::plots for valid names.", 
         call. = FALSE)
  
  return(fun)
}
