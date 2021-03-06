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

.printfr <- function(x, digits, ...) {
  print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
}
.median_and_madsd <- function(x) {
  cbind(Median = apply(x, 2, median), MAD_SD = apply(x, 2, mad))
}

#' Print method for stanjm objects
#' 
#' The \code{print} method for \code{\link[=stanjm-object]{stanjm}}
#' objects displays a compact summary of a joint model fitted using the 
#' \code{\link{stan_jm}} modelling function. For more detailed summary statistics 
#' and diagnostics use the \code{\link[=summary.stanjm]{summary}} method.
#' 
#' @export
#' @method print stanjm
#' @templateVar stanjmArg x
#' @template args-stanjm-object
#' @param digits Number of digits to use for formatting numbers.
#' @param ... Ignored.
#' @return Returns \code{x}, invisibly.
#' @details 
#' The output from the \code{print} method includes very brief summary statistics 
#' for the fixed effects parameter estimates (described below, and which are 
#' displayed separately for each longitudinal or event submodel), the additional
#' parameters such as the residual error standard deviation, dispersion parameter, 
#' or the Weibull shape parameter or baseline hazard coefficients. Also displayed 
#' are the standard deviation and correlation parameter estimates for the random effects.  
#'
#' \subsection{Point estimates}{
#' Point estimates are medians computed from the posterior sample generated 
#' via the MCMC simulations. The point estimates reported are the same as the 
#' values returned by the \code{\link[=coef.stanjm]{coef}} method.
#' }
#' \subsection{Uncertainty estimates}{
#' The standard deviations (labeled MAD_SD in the print output) are computed from
#' the same posterior sample used to calculate the point estimates and are proportional 
#' to the median absolute deviation (\code{\link[stats]{mad}}) from the median. 
#' Compared to the raw posterior standard deviation, the MAD_SD will be more 
#' robust for long-tailed distributions. These are the same as the values 
#' returned by \code{\link[=se.stanjm]{se}}.
#' } 
#' 
#' @seealso \code{\link{summary.stanjm}}, \code{\link{stanjm-methods}}
#'
print.stanjm <- function(x, digits = 3, ...) {
  print(x$call) 
  
  M <- x$n_markers
  link    <- sapply(1:M, function(m) x$family[[m]]$link)

  mat <- as.matrix(x$stanfit)
  nms <- collect_nms(rownames(x$stan_summary), M, value = TRUE)
      
  # Estimates table for longitudinal submodel(s)
  for (m in 1:M) {
    terms_m <- terms(x)[[m]]
    sel <- attr(terms_m, "response")
    yvar <- rownames(attr(terms_m, "factors"))[sel]
    cat(paste0("\nLongitudinal submodel", if (M > 1) paste0(" ", m), 
               ": ", yvar,"\n"))
    coef_mat <- mat[, c(nms$y[[m]], nms$y_extra[[m]]), drop = FALSE]
    
    # Calculate median and MAD
    estimates <- .median_and_madsd(coef_mat)
    
    # Add column with eform
    if (link[m] %in% c("log", "logit")) 
      estimates <- cbind(estimates, 
        "exp(Median)" = c(exp(estimates[nms$y[[m]], "Median"]), 
                          rep(NA, length(nms$y_extra[[m]]))))
    
    # Print estimates
    rownames(estimates) <- 
      gsub(paste0("^Long", m, "\\|"), "", rownames(estimates))     
    .printfr(estimates, digits, ...)
  }
  
  # Estimates table for event submodel
    cat("\nEvent submodel:\n")   
    coef_mat <- mat[, c(nms$e, nms$a, nms$e_extra), drop = FALSE]
    
    # Calculate median and MAD
    estimates <- .median_and_madsd(coef_mat)
  
    # Add column with eform
    estimates <- cbind(estimates, 
      "exp(Median)" = c(exp(estimates[c(nms$e, nms$a), "Median"]), 
                        rep(NA, length(nms$e_extra))))
  
    rownames(estimates) <- gsub("^Event\\|", "", rownames(estimates))  
    rownames(estimates) <- gsub("^Assoc\\|", "", rownames(estimates))   
    .printfr(estimates, digits, ...)

  # Estimates table for group-level random effects
  cat("\nGroup-level random effects:\n") 
  print(VarCorr(x), digits = digits + 1, ...)
  cat("Num. levels:", paste(names(ngrps(x)), unname(ngrps(x)), 
                            collapse = ", "), "\n")  
  
  invisible(x)
}


#' Summary method for stanjm objects
#' 
#' The \code{summary} method for \code{\link[=stanjm-object]{stanjm}}
#' objects returns a summary of parameter estimates and MCMC convergence diagnostics 
#' (Monte Carlo error, effective sample size, Rhat) for a joint model fitted
#' using the \code{\link{stan_jm}} modelling function.
#'  
#' @export
#' @method summary stanjm
#' 
#' @templateVar stanjmArg object
#' @template args-stanjm-object
#' @template args-regex-pars
#' 
#' @param ... Currently ignored.
#' @param pars An optional character vector specifying a subset of parameters to
#'   display. Parameters can be specified by name or several shortcuts can be 
#'   used. 
#'   Using \code{pars = "long"} will display the 
#'   parameter estimates for the longitudinal submodels only (excluding random
#'   effects, but including dispersion parameters).
#'   Using \code{pars = "event"} will display the 
#'   parameter estimates for the event submodel only, including any association
#'   parameters. 
#'   Using \code{pars = "assoc"} will display only the 
#'   association parameters. 
#'   Using \code{pars = "fixef"} will display all fixed effects, but not
#'   the random effects or the additional parameters such as dispersion, etc. 
#'   Using \code{pars = "beta"} will display only the 
#'   fixed effect regression coefficients (excluding the intercept terms).
#'   Using \code{pars = "alpha"} will display only the 
#'   fixed effect intercept terms.
#'   The estimates for the random effects can be selected using \code{pars = "b"}
#'   or \code{pars = "varying"}.
#'   See \strong{Examples}.
#'   If both \code{pars} and \code{regex_pars} are set to \code{NULL} then all 
#'   fixed effect regression coefficients are selected, as well as any 
#'   additional parameters such as the residual error, dispersion,
#'   or baseline hazard parameters, and the log posterior. 
#' @param probs An optional numeric vector of probabilities passed to 
#'   \code{\link[stats]{quantile}}, for calculating the quantiles of the 
#'   posterior distribution for each parameter.
#' @param digits Number of digits to use for formatting numbers when printing. 
#'   When calling \code{summary}, the value of digits is stored as the 
#'   \code{"print.digits"} attribute of the returned object.
#'   
#' @return The \code{summary} method returns an object of class 
#'   \code{"summary.stanjm"}, which is a matrix of summary statistics and 
#'   diagnostics, with additional information stored as attributes.
#'   The \code{print} method for \code{summary.stanjm}
#'   objects is called for its side effect and just returns its input. The 
#'   \code{as.data.frame} method for \code{summary.stanjm} objects converts 
#'   the matrix to a data.frame, preserving row and column names but dropping  
#'   the \code{print}-related attributes.
#' 
#' @seealso \code{\link{print.stanjm}}, \code{\link{stanjm-methods}}
#' 
#' @examples
#' if (!exists("examplejm")) example(examplejm)
#' 
#' # Only showing 10th and 90th percentile   
#' summary(examplejm, probs = c(0.1, 0.9))
#' 
#' # These produce the same output for this example, 
#' # but the first method can be used for any model
#' summary(examplejm, pars = c("long"))
#' summary(examplejm, pars = c("Long1|(Intercept)",
#'                              "Long1|year", 
#'                              "Long1|sigma")) 
#' 
#' # Only show parameters for event submodel
#' summary(examplejm, pars = "event")
#'
#' # Only show the association parameter for the current
#' # value of the linear predictor from the longitudinal submodel
#' summary(examplejm, pars = "Assoc|Long1:eta-value")
#' # or since there is only one association parameter in the 
#' # model we can just use the following shortcut
#' summary(examplejm, pars = "assoc")
#'   
#' # Only show random effects parameters
#' summary(examplejm, pars = "b")
#' as.data.frame(summary(examplejm, pars = "b"))
#'  
#' # To obtain only the random intercepts we could also use a
#' # regular expression in the regex argument as follows
#' summary(examplejm, pars = NULL,
#'         regex_pars = "\\(Intercept\\)\\sid")
#'  
#' @importMethodsFrom rstan summary
#' 
summary.stanjm <- function(object, pars = NULL, regex_pars = NULL, 
                            probs = NULL, digits = 3, ...) {
  pars <- collect_pars(object, pars, regex_pars)
  M <- object$n_markers
 
  # Outcome variable for each submodel
  y_vars <- sapply(1:M, function(m, object) {
    terms_m <- terms(object)[[m]]
    sel <- attr(terms_m, "response")
    ret <- rownames(attr(terms_m, "factors"))[sel]
  }, object = object)
   
  # Family and link for each submodel
  fam <- sapply(object$family, function(x) 
           paste0(x$family, " (", x$link, ")")) 
  
  # Construct summary table  
  args <- list(object = object$stanfit)
  if (!is.null(probs)) 
    args$probs <- probs
  out <- do.call("summary", args)$summary
  
  nms <- collect_nms(rownames(object$stan_summary), M, value = TRUE)
  if (!is.null(pars)) {
    pars2 <- NA     
    if ("alpha" %in% pars) pars2 <- c(pars2, nms$alpha)
    if ("beta" %in% pars) pars2 <- c(pars2, nms$beta)
    if ("long" %in% pars) pars2 <- c(pars2, unlist(nms$y), unlist(nms$y_extra))
    if ("event" %in% pars) pars2 <- c(pars2, nms$e, nms$a, nms$e_extra)
    if ("assoc" %in% pars) pars2 <- c(pars2, nms$a)      
    if ("fixef" %in% pars) pars2 <- c(pars2, unlist(nms$y), nms$e, nms$a)
    if ("b" %in% pars) pars2 <- c(pars2, nms$b)
    pars2 <- c(pars2, setdiff(pars, 
                              c("alpha", "beta", "varying", "b",
                                "long", "event", "assoc", "fixef")))
    pars <- pars2[!is.na(pars2)]
  } else {
    pars <- c(unlist(nms$y), unlist(nms$y_extra),
              nms$e, nms$a, nms$e_extra, "log-posterior")
  }
  
  out <- out[rownames(out) %in% pars, , drop = FALSE]
  out <- out[!grepl(":_NEW_", rownames(out), fixed = TRUE), , drop = FALSE]
  stats <- colnames(out)
  if ("n_eff" %in% stats)
    out[, "n_eff"] <- round(out[, "n_eff"])
  if ("se_mean" %in% stats) # So people don't confuse se_mean and sd
    colnames(out)[stats %in% "se_mean"] <- "mcse"
    
  # Reorder rows of output table
  nms_tmp <- rownames(out)  
  nms_tmp_y <- lapply(1:M, function(m) 
                 grep(paste0("^Long", m, "\\|"), nms_tmp, value = TRUE))
  nms_tmp_e <- grep("^Event\\|", nms_tmp, value = TRUE)
  nms_tmp_a <- grep("^Assoc\\|", nms_tmp, value = TRUE)
  nms_tmp_b <- b_names(nms_tmp, value = TRUE)
  nms_tmp_lp <- grep("^log-posterior$", nms_tmp, value = TRUE)
  out <- out[c(unlist(nms_tmp_y), nms_tmp_e, nms_tmp_a, nms_tmp_b, nms_tmp_lp), , drop = FALSE]
  
  # Run times
  times <- round((rstan::get_elapsed_time(object$stanfit))/60, digits = 1)
  times <- cbind(times, total = rowSums(times))

  # Output object
  structure(out, 
            call = object$call, 
            algorithm = object$algorithm,
            n_markers = object$n_markers,
            n_subjects = object$n_subjects,
            n_grps = object$n_grps,
            n_events = object$n_events,
            n_yobs = object$n_yobs,
            y_vars = y_vars,
            id_var = object$id_var,
            time_var = object$time_var,
            family = fam,
            base_haz = object$base_haz$type,
            posterior_sample_size = posterior_sample_size(object),
            times = times,
            print.digits = digits, 
            class = "summary.stanjm")
}

#' @rdname summary.stanjm
#' @export
#' @method print summary.stanjm
#'
#' @param x An object of class \code{"summary.stanjm"}.
#'
#' @importMethodsFrom rstan get_elapsed_time
print.summary.stanjm <- function(x, digits = max(1, attr(x, "print.digits")), 
                                  ...) {
  atts <- attributes(x)
  M <- atts$n_markers
  
  print(atts$call)
  
  cat(paste0("\n", if (M == 1) "Uni" else "Multi", 
             "variate joint model, consisting of:")) 
  for (m in 1:M) {
    cat(paste0("\n  Outcome variable", 
               if (M > 1) paste0(" (Long ", m, ")"), 
               ": ", atts$y_vars[m]))
  }
  for (m in 1:M) {
    cat(paste0("\n  Family", 
               if (M > 1) paste0(" (Long ", m, ")"), 
               ": ", atts$family[m]))
  }
  cat(paste0("\n  Baseline hazard: ", atts$base_haz))  
  cat(paste0("\n  Clustering variables: ", paste(names(atts$n_grps), sep = ",")))
  if (!is.null(atts$n_subjects))
    cat(paste0("\n  Num. subjects (", atts$id_var, "): ", atts$n_subjects))
  cat(paste0("\n  Num. events: ", atts$n_events, " (", 
            round(100 * atts$n_events/atts$n_subjects, 1), "%)"))
  cat("\n  Num. long observations: ")
  cat(paste0(atts$n_yobs, if (M > 1) paste0(" (Long ", 1:M, ")"), collapse = ", "))
  cat("\n  Posterior sample size:", atts$posterior_sample_size, "MCMC iterations")
  
  cat("\n\nTime taken for sampling (mins):\n")
  print(atts$times)
  
  cat("\nEstimates:\n")
  sel <- which(colnames(x) %in% c("mcse", "n_eff", "Rhat"))
  if (!length(sel)) {
    .printfr(x, digits)
  } else {
    xtemp <- x[, -sel, drop = FALSE]
    colnames(xtemp) <- paste(" ", colnames(xtemp))
    .printfr(xtemp, digits)
    cat("\nDiagnostics:\n")
    mcse_rhat <- format(round(x[, c("mcse", "Rhat"), drop = FALSE], digits), 
                        nsmall = digits)
    n_eff <- format(x[, "n_eff", drop = FALSE], drop0trailing = TRUE)
    print(cbind(mcse_rhat, n_eff), quote = FALSE)
    cat("\nFor each parameter, mcse is Monte Carlo standard error, ", 
        "n_eff is a crude measure of effective sample size, ", 
        "and Rhat is the potential scale reduction factor on split chains", 
        " (at convergence Rhat=1).\n", sep = '')
  }
  invisible(x)
}

#' @rdname summary.stanjm
#' @method as.data.frame summary.stanjm
#' @export
as.data.frame.summary.stanjm <- function(x, ...) {
  as.data.frame(unclass(x), ...)
}
