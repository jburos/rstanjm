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

# Issue warning if high rhat values
# 
# @param rhats Vector of rhat values.
# @param threshold Threshold value. If any rhat values are above threshold a 
#   warning is issued.
check_rhats <- function(rhats, threshold = 1.1) {
  if (any(rhats > threshold, na.rm = TRUE)) 
    warning("Markov chains did not converge! Do not analyze results!", 
            call. = FALSE, noBreaks. = TRUE)
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

# If a is NULL (and Inf, respectively) return b, otherwise just return a
# @param a,b Objects
`%ORifNULL%` <- function(a, b) {
  if (is.null(a)) b else a
}
`%ORifINF%` <- function(a, b) {
  if (a == Inf) b else a
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

# Convert 2-level factor to 0/1
fac2bin <- function(y) {
  if (!is.factor(y)) 
    stop("Bug found: non-factor as input to fac2bin.", 
         call. = FALSE)
  if (!identical(nlevels(y), 2L)) 
    stop("Bug found: factor with nlevels != 2 as input to fac2bin.", 
         call. = FALSE)
  as.integer(y != levels(y)[1L])
}

# Test for a given family
#
# @param x A character vector (probably x = family(fit)$family)
is.binomial <- function(x) x == "binomial"
is.gaussian <- function(x) x == "gaussian"
is.gamma <- function(x) x == "Gamma"
is.ig <- function(x) x == "inverse.gaussian"
is.nb <- function(x) x == "neg_binomial_2"
is.poisson <- function(x) x == "poisson"

is_scobit <- function(object) {
  validate_stanreg_object(object)
  if (!is(object, "polr")) return(FALSE)
  return("alpha" %in% rownames(object$stan_summary))
}

# Check family argument
#
# @param f The \code{family} argument specified by user (or the default).
# @return If no error is thrown, then either \code{f} itself is returned (if
#   already a family) or the family object created from \code{f} is returned (if
#   \code{f} is a string or function).
validate_family <- function(f) {
  if (is.character(f)) 
    f <- get(f, mode = "function", envir = parent.frame(2))
  if (is.function(f)) 
    f <- f()
  if (!is(f, "family")) 
    stop("'family' must be a family.", call. = FALSE)
  
  return(f)
}

# Make prior.info list
# @param user_call The user's call, i.e. match.call(expand.dots = TRUE).
# @param function_formals Formal arguments of the stan_* function, i.e.
#   formals().
# @return A list containing information about the prior distributions and
#   options used.
get_prior_info <- function(user_call, function_formals) {
  user <- grep("prior", names(user_call), value = TRUE)
  default <- setdiff(grep("prior", names(function_formals), value = TRUE), 
                     user)
  U <- length(user)
  D <- length(default)
  priors <- list()
  for (j in 1:(U + D)) {
    if (j <= U) {
      priors[[user[j]]] <- try(eval(user_call[[user[j]]]), silent = TRUE)
    } else {
      priors[[default[j-U]]] <- try(eval(function_formals[[default[j-U]]]), 
                                    silent = TRUE)
    } 
  }
  
  return(priors)
}

# Return names of the last dimension in a matrix/array (e.g. colnames if matrix)
#
# @param x A matrix or array
last_dimnames <- function(x) {
  ndim <- length(dim(x))
  dimnames(x)[[ndim]]
}

# Methods for creating linear predictor
#
# Make linear predictor vector from x and point estimates for beta, or linear 
# predictor matrix from x and full posterior sample of beta.
#
# @param beta A vector or matrix or parameter estimates.
# @param x Predictor matrix.
# @param offset Optional offset vector.
# @return A vector or matrix.
linear_predictor <- function(beta, x, offset = NULL) {
  UseMethod("linear_predictor")
}
linear_predictor.default <- function(beta, x, offset = NULL) {
  eta <- as.vector(if (NCOL(x) == 1L) x * beta else x %*% beta)
  if (length(offset))
    eta <- eta + offset
  
  return(eta)
}
linear_predictor.matrix <- function(beta, x, offset = NULL) {
  if (NCOL(beta) == 1L) 
    beta <- as.matrix(beta)
  eta <- beta %*% t(x)
  if (length(offset)) 
    eta <- sweep(eta, 2L, offset, `+`)
  
  return(eta)
}

# Wrapper for rstan::summary
# @param stanfit A stanfit object created using rstan::sampling or rstan::vb
# @return A matrix of summary stats
make_stan_summary <- function(stanfit) {
  levs <- c(0.5, 0.8, 0.95, 0.99)
  qq <- (1 - levs) / 2
  probs <- sort(c(0.5, qq, 1 - qq))
  rstan::summary(stanfit, probs = probs, digits = 10)$summary  
}

# Get the correct column name to use for selecting the median
#
# @param algorithm String naming the estimation algorithm (probably
#   \code{fit$algorithm}).
# @return Either \code{"50%"} or \code{"Median"} depending on \code{algorithm}.
select_median <- function(algorithm) {
  switch(algorithm, 
         sampling = "50%",
         meanfield = "50%",
         fullrank = "50%",
         optimizing = "Median",
         stop("Bug found (incorrect algorithm name passed to select_median)", 
              call. = FALSE))
}

# Maybe broadcast 
#
# @param x A vector or scalar.
# @param n Number of replications to possibly make. 
# @return If \code{x} has no length the \code{0} replicated \code{n} times is
#   returned. If \code{x} has length 1, the \code{x} replicated \code{n} times
#   is returned. Otherwise \code{x} itself is returned.
maybe_broadcast <- function(x, n) {
  if (!length(x)) {
    rep(0, times = n)
  } else if (length(x) == 1L) {
    rep(x, times = n)
  } else {
    x
  }
}

# Create a named list using specified names or, if names are omitted, using the
# names of the objects in the list
#
# @param ... Objects to include in the list. 
# @return A named list.
nlist <- function(...) {
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names) FALSE else nzchar(names(out))
  if (all(has_name)) 
    return(out)
  nms <- as.character(m)[-1L]
  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  } 
  
  return(out)
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

# Check and set scale parameters for priors
#
# @param scale Value of scale parameter (can be NULL).
# @param default Default value to use if \code{scale} is NULL.
# @param link String naming the link function.
# @return If a probit link is being used, \code{scale} (or \code{default} if
#   \code{scale} is NULL) is scaled by \code{dnorm(0) / dlogis(0)}. Otherwise
#   either \code{scale} or \code{default} is returned.
set_prior_scale <- function(scale, default, link) {
  stopifnot(is.numeric(default), is.character(link))
  if (is.null(scale)) 
    scale <- default
  if (link == "probit")
    scale <- scale * dnorm(0) / dlogis(0)
  
  return(scale)
}

# Check for positive scale or df parameter (NULL ok)
#
# @param x The value to check.
# @return Either an error is thrown or \code{TRUE} is returned invisibly.
validate_parameter_value <- function(x) {
  nm <- deparse(substitute(x))
  if (!is.null(x)) {
    if (!is.numeric(x)) 
      stop(nm, " should be NULL or numeric", call. = FALSE)
    if (any(x <= 0)) 
      stop(nm, " should be positive", call. = FALSE)
  }
  invisible(TRUE)
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


#----------  Functions from: stanreg-methods.R ----------# 

justRE <- function(f, response = FALSE) {
  response <- if (response && length(f) == 3) f[[2]] else NULL
  reformulate(paste0("(", vapply(lme4::findbars(f), 
                                 function(x) paste(deparse(x, 500L), 
                                                   collapse = " "), 
                                 ""), ")"), 
              response = response)
}


#----------  Functions from: as.matrix.stanreg.R ----------# 

STOP_no_draws <- function() stop("No draws found.", call. = FALSE)

check_missing_pars <- function(x, pars) {
  notfound <- which(!pars %in% last_dimnames(x))
  if (length(notfound)) 
    stop(
      "No parameter(s) ", 
      paste(pars[notfound], collapse = ", "), 
      call. = FALSE
    )
}

exclude_lp_and_ppd <- function(pars) {
  grep(
    pattern = "mean_PPD|log-posterior", 
    x = pars, 
    invert = TRUE, 
    value = TRUE
  )
}


#----------  Functions from: posterior_predict.R ----------# 

pp_b_ord <- function(b, Z_names) {
  ord <- sapply(Z_names, FUN = function(x) {
    m <- grep(paste0("b[", x, "]"), colnames(b), fixed = TRUE)
    len <- length(m)
    if (len == 1) 
      return(m)
    if (len > 1) 
      stop("multiple matches bug")
    m <- grep(paste0("b[", sub(" (.*):.*$", " \\1:_NEW_\\1", x), "]"), 
              colnames(b), fixed = TRUE)
    if (len == 1)
      return(m)
    if (len > 1)
      stop("multiple matches bug")
    x <- strsplit(x, split = ":", fixed = TRUE)[[1]]
    stem <- strsplit(x[[1]], split = " ", fixed = TRUE)[[1]]
    x <- paste(x[1], x[2], paste0("_NEW_", stem[2]), x[2], sep = ":")
    m <- grep(paste0("b[", x, "]"), colnames(b), fixed = TRUE)
    len <- length(m)
    if (len == 1)
      return(m)
    if (len > 1)
      stop("multiple matches bug")
    x <- paste(paste(stem[1], stem[2]), paste0("_NEW_", stem[2]), sep = ":")
    m <- grep(paste0("b[", x, "]"), colnames(b), fixed = TRUE)
    len <- length(m)
    if (len == 1)
      return(m)
    if (len > 1)
      stop("multiple matches bug")
    stop("no matches bug")    
  })
  b[, ord, drop = FALSE]
}


#----------  Functions from: pp_check.R ----------# 

#' @importFrom ggplot2 geom_hline geom_point geom_path labs facet_wrap ggtitle
pp_check_binned_resid <- function(object, n = 1, ...) {
  if (!requireNamespace("arm", quietly = TRUE)) 
    stop("This plot requires the 'arm' package (install.packages('arm'))", 
         call. = FALSE)
  
  binner <- function(rep_id, ey, r, nbins) {
    br <- arm::binned.resids(ey, r, nbins)$binned[, c("xbar", "ybar", "2se")]
    if (length(dim(br)) < 2L)
      br <- t(br)
    colnames(br) <- c("xbar", "ybar", "se2")
    data.frame(rep = paste0("yrep_", rep_id), br)
  }
  
  dat <- pp_data(object, newdata = NULL)
  stanmat <- as.matrix.stanreg(object)
  sel <- sample(nrow(stanmat), size = n)
  beta <- stanmat[sel, 1:ncol(dat$x), drop = FALSE]
  eta <- linear_predictor(beta, dat$x, dat$offset)
  inverse_link <- linkinv(object)
  Ey <- inverse_link(eta)
  y <- get_y(object)
  if (NCOL(y) == 2L) 
    y <- y[, 1L] / rowSums(y)
  ytmp <- if (is.factor(y)) 
    fac2bin(y) else y
  resids <- sweep(-Ey, MARGIN = 2L, STATS = ytmp, "+")
  ny <- length(y)
  stopifnot(ny == ncol(Ey))
  if (ny >= 100) {
    nbins <- floor(sqrt(ny))
  } else if (ny > 10 && ny < 100) {
    nbins <- 10
  } else { # nocov start
    # if (ny <= 10)
    nbins <- floor(ny / 2) 
  } # nocov end
  
  binned <- binner(rep_id = 1, ey = Ey[1L, ], r = resids[1L, ], nbins)
  if (n > 1) {
    for (i in 2:nrow(resids))
      binned <- rbind(binned, binner(rep_id = i, ey = Ey[i, ], 
                                     r = resids[i, ], nbins))
  }
  dots <- list(...)
  line_color <- dots$color %ORifNULL% .PP_FILL
  line_size <- dots$size %ORifNULL% 1
  pt_color <- dots$fill %ORifNULL% .PP_VLINE_CLR
  base <- ggplot(binned, aes_string(x = "xbar"))
  graph <- base + 
    geom_hline(yintercept = 0, linetype = 2) + 
    geom_path(aes_string(y = "se2"), color = line_color, size = line_size) + 
    geom_path(aes_string(y = "-se2"), color = line_color, size = line_size) + 
    geom_point(aes_string(y = "ybar"), shape = 19, color = pt_color) + 
    labs(x = "Expected Values", y = "Average Residual \n (with 2SE bounds)") + 
    ggtitle("Binned Residuals")
  
  if (n > 1) 
    graph <- graph + facet_wrap(~rep, scales = "free")
  
  return(graph + pp_check_theme(no_y = FALSE))
}

as.matrix.stanreg <- function(x, ..., pars = NULL, regex_pars = NULL) {
  pars <- collect_pars(x, pars, regex_pars)
  user_pars <- !is.null(pars)
  
  if (used.optimizing(x)) {
    mat <- x$asymptotic_sampling_dist
    if (is.null(mat)) 
      STOP_no_draws()
    if (!user_pars) {
      dispersion <- c("sigma", "scale", "shape", "lambda", "overdispersion")
      pars <- c(names(coef(x)), # return with coefficients first
                dispersion[which(dispersion %in% colnames(mat))])
    }
  } else { # used mcmc or vb
    mat <- as.matrix(x$stanfit)
    if (!user_pars)
      pars <- exclude_lp_and_ppd(colnames(mat))
  }
  if (user_pars)
    check_missing_pars(mat, pars)
  
  mat <- mat[, pars, drop = FALSE]
  if (!is.mer(x))
    return(mat)
  unpad_reTrms(mat)
}

set_geom_args <- function(defaults, ...) {
  dots <- list(...)
  if (!length(dots)) 
    return(defaults)
  dot_names <- names(dots)
  def_names <- names(defaults)
  for (j in seq_along(def_names)) {
    if (def_names[j] %in% dot_names)
      defaults[[j]] <- dots[[def_names[j]]]
  }
  extras <- setdiff(dot_names, def_names)
  if (length(extras)) {
    for (j in seq_along(extras))
      defaults[[extras[j]]] <- dots[[extras[j]]]
  }
  return(defaults)
}

#' @importFrom ggplot2 ggtitle element_blank element_line element_text
#'   theme_classic theme %+replace%
pp_check_theme <- function(no_y = TRUE) {
  blank <- element_blank()
  thm <- theme_classic() +
    theme(axis.line = element_line(color = "#222222"),
          axis.line.y = if (no_y) blank  else element_line(size = 0.5),
          axis.line.x = element_line(size = 2),
          axis.title = element_text(face = "bold", size = 13),
          strip.background = element_blank(),
          strip.text = element_text(color = "black", face = "bold"),
          legend.position = "none",
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 13),
          plot.title = element_text(size = 18))
  if (no_y)
    thm <- thm %+replace% theme(axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank())
  return(thm)
}
