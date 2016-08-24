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

# Create a stanjm object
#
# @param object A list provided by one of the \code{stan_*} modeling functions.
# @return A stanjm object.
#
stanjm <- function(object) {
  opt        <- object$algorithm == "optimizing"
  mer        <- rep(1L, object$M)
  stanfit    <- object$stanfit
  M          <- object$M
  family     <- object$family
  assoc      <- object$assoc
  y          <- object$y
  x          <- object$x
  xq         <- object$xq
  dxdtq      <- object$dxdtq  
  e_x        <- object$e_x
  eventtime  <- object$eventtime
  d          <- object$d  
  dimensions <- object$dimensions
  y_cnms <- object$y_cnms
  y_flist <- object$y_flist
  y_nms      <- lapply(y, names)
  if (opt) {
    stop("Optimisation not implemented for stan_jm")
  } else {
    stan_summary <- rstanarm:::make_stan_summary(stanfit)
    nms <- collect_nms(rownames(stan_summary), M)
    
    # Coefs and SEs for longitudinal submodel(s)                    
    y_coefs <- lapply(1:M, function(m)
      stan_summary[c(nms$y[[m]], nms$y_b[[m]]), rstanarm:::select_median(object$algorithm)])
    y_stanmat <- lapply(1:M, function(m) 
      as.matrix(stanfit)[, c(nms$y[[m]], nms$y_b[[m]]), drop = FALSE])
    y_ses <- lapply(y_stanmat, function(m) apply(m, 2L, mad))
    y_covmat <- lapply(y_stanmat, cov)
    for (m in 1:M) {
      rownames(y_covmat[[m]]) <- colnames(y_covmat[[m]]) <- rownames(stan_summary)[c(nms$y[[m]], nms$y_b[[m]])]
    }
 
    # Coefs and SEs for event submodel    
    e_coefs <- stan_summary[c(nms$e, nms$a), rstanarm:::select_median(object$algorithm)]        
    if (length(e_coefs) == 1L) names(e_coefs) <- rownames(stan_summary)[c(nms$e, nms$a)[1L]]
    e_stanmat <- as.matrix(stanfit)[, c(nms$e, nms$a), drop = FALSE]
    e_ses <- apply(e_stanmat, 2L, mad)    
    e_covmat <- cov(e_stanmat)
    rownames(e_covmat) <- colnames(e_covmat) <- rownames(stan_summary)[c(nms$e, nms$a)]

    if (object$algorithm == "sampling") 
      rstanarm:::check_rhats(stan_summary[, "Rhat"])
  }
  
  # Linear predictor, fitted values
  y_eta <- lapply(1:M, function(m) rstanarm:::linear_predictor.default(y_coefs[[m]], x[[m]], object$offset))
  y_mu  <- lapply(1:M, function(m) family[[m]]$linkinv(y_eta[[m]]))

  # Residuals
  y_tmp <- lapply(1:M, function(m) if (is.factor(y[[m]])) rstanarm:::fac2bin(y[[m]]) else y[[m]])
  y_residuals <- lapply(1:M, function(m) y_tmp[[m]] - y_mu[[m]])
  for (m in 1:M) {
    names(y_eta[[m]]) <- names(y_mu[[m]]) <- names(y_residuals[[m]]) <- y_nms[[m]]
  }

  # Remove padding
  y_coefs <- lapply(y_coefs, unpad_reTrms.default)
  y_ses   <- lapply(y_ses, unpad_reTrms.default)

  out <- rstanarm:::nlist(
    coefficients = list_nms(c(y_coefs, list(e_coefs)), M), 
    ses = list_nms(c(y_ses, list(e_ses)), M),
    fitted.values = list_nms(y_mu, M),
    linear.predictors = list_nms(y_eta, M),
    residuals = list_nms(y_residuals, M), 
    df.residual = if (opt) df.residual else NA_integer_, 
    covmat = list_nms(c(y_covmat, list(e_covmat)), M),
    n_markers = M,
    n_subjects = object$Npat,
    n_grps = object$n_grps,
    n_events = sum(d > 0),
    n_yobs = object$y_N,
    id_var = object$id_var,
    time_var = object$time_var,
    cnms = object$cnms, 
    y_cnms = list_nms(y_cnms, M), 
    y_flist = list_nms(y_flist, M),
    assoc,
    y = list_nms(y, M), 
    x = list_nms(x, M), 
    xq = list_nms(xq, M), 
    dxdtq = if (sum(assoc$etaslope)) list_nms(dxdtq, M) else NULL,
    eventtime, d,     
    model = object$model, 
    dataLong = object$dataLong,
    dataEvent = object$dataEvent,     
    family = list_nms(family, M), 
    base_haz = object$base_haz,
#    offset = if (any(object$offset != 0)) object$offset else NULL,
#    weights = object$weights, 
#    prior.weights = object$weights, 
#    contrasts = object$contrasts, 
    na.action = object$na.action,
    call = object$call, 
    formula = list_nms(object$formula, M), 
    terms = object$terms,
    prior.info = object$prior.info,
    algorithm = object$algorithm,
    stan_summary,  
    stanfit = if (opt) stanfit$stanfit else stanfit,
    glmod = object$glmod
  )
  if (opt) 
    out$asymptotic_sampling_dist <- stanmat
  
  structure(out, class = c("stanjm", "stanreg", "lmerMod"))
}
