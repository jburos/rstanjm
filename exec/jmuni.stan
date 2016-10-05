/*
    This file is part of rstanarm.
    Copyright (C) 2015, 2016 Trustees of Columbia University
    
    rstanarm is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    rstanarm is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with rstanarm.  If not, see <http://www.gnu.org/licenses/>.
*/

functions {

  /** 
   * Apply inverse link function to linear predictor
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_gauss(vector eta, int link) {
    if (link < 1 || link > 3) reject("Invalid link");
    if (link < 3)  # link = identity or log 
      return eta; # return eta for log link too bc will use lognormal
    else {# link = inverse
      vector[rows(eta)] mu;
      for(n in 1:rows(eta)) mu[n] = inv(eta[n]); 
      return mu;
    }
  }
  
  /** 
  * Apply inverse link function to linear predictor
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_gamma(vector eta, int link) {
    if (link < 1 || link > 3) reject("Invalid link");
    if (link == 1)  return eta;
    else if (link == 2) return exp(eta);
    else {
      vector[rows(eta)] mu;
      for(n in 1:rows(eta)) mu[n] = inv(eta[n]); 
      return mu;
    }
  }
  
  /** 
  * Apply inverse link function to linear predictor
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_inv_gaussian(vector eta, int link) {
    if (link < 1 || link > 4) reject("Invalid link");
    if (link == 1)  return eta;
    else if (link == 2) return exp(eta);
    else {
      vector[rows(eta)] mu;
      if (link == 3) for( n in 1:rows(eta)) mu[n] = inv(eta[n]);
      else for (n in 1:rows(eta)) mu[n] = inv_sqrt(eta[n]);      
      return mu;
    }
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_gauss(vector y, vector eta, real sigma, int link) {
    vector[rows(eta)] ll;
    if (link < 1 || link > 3) reject("Invalid link");
    if (link == 2) # link = log
      for (n in 1:rows(eta)) ll[n] = lognormal_lpdf(y[n] | eta[n], sigma);
    else { # link = idenity or inverse
      vector[rows(eta)] mu;
      mu = linkinv_gauss(eta, link);
      for (n in 1:rows(eta)) ll[n] = normal_lpdf(y[n] | mu[n], sigma);
    }
    return ll;
  }
  
  real GammaReg(vector y, vector eta, real shape, 
                int link, real sum_log_y) {
    real ret;
    if (link < 1 || link > 3) reject("Invalid link");
    ret = rows(y) * (shape * log(shape) - lgamma(shape)) +
      (shape - 1) * sum_log_y;
    if (link == 2)      # link is log
      ret = ret - shape * sum(eta) - shape * sum(y ./ exp(eta));
    else if (link == 1) # link is identity
      ret = ret - shape * sum(log(eta)) - shape * sum(y ./ eta);
    else                # link is inverse
      ret = ret + shape * sum(log(eta)) - shape * dot_product(eta, y);
    return ret;
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_gamma(vector y, vector eta, real shape, int link) {
    vector[rows(eta)] ll;
    if (link < 1 || link > 3) reject("Invalid link");
    if (link == 3) { # link = inverse
      for (n in 1:rows(eta)) {
        ll[n] = gamma_lpdf(y[n] | shape, shape * eta[n]);
      }
    }
    else if (link == 2) { # link = log
      for (n in 1:rows(eta)) {
        ll[n] = gamma_lpdf(y[n] | shape, shape / exp(eta[n]));
      }
    }
    else { # link = identity
      for (n in 1:rows(eta)) {
        ll[n] = gamma_lpdf(y[n] | shape, shape / eta[n]);
      }
    }
    return ll;
  }
  
  /** 
  * inverse Gaussian log-PDF (for data only, excludes constants)
  *
  * @param y The vector of outcomes
  * @param eta The vector of linear predictors
  * @param lambda A positive scalar nuisance parameter
  * @param link An integer indicating the link function
  * @return A scalar
  */
  real inv_gaussian(vector y, vector mu, real lambda, 
                    real sum_log_y, vector sqrt_y) {
    return 0.5 * rows(y) * log(lambda / (2 * pi())) - 
      1.5 * sum_log_y - 
      0.5 * lambda * dot_self( (y - mu) ./ (mu .* sqrt_y) );
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param eta The linear predictors
  * @param lamba A positive scalar nuisance parameter
  * @param link An integer indicating the link function
  * @param log_y A precalculated vector of the log of y
  * @param sqrt_y A precalculated vector of the square root of y
  * @return A vector of log-likelihoods
  */
  vector pw_inv_gaussian(vector y, vector eta, real lambda, 
                         int link, vector log_y, vector sqrt_y) {
    vector[rows(y)] ll;
    vector[rows(y)] mu;
    if (link < 1 || link > 4) reject("Invalid link");
    mu = linkinv_inv_gaussian(eta, link);
    for (n in 1:rows(y))
      ll[n] = -0.5 * lambda * square( (y[n] - mu[n]) / (mu[n] * sqrt_y[n]) );
    ll = ll + 0.5 * log(lambda / (2 * pi())) - 1.5 * log_y;
    return ll;
  }
  
  /** 
  * PRNG for the inverse Gaussian distribution
  *
  * Algorithm from wikipedia 
  *
  * @param mu The expectation
  * @param lambda The dispersion
  * @return A draw from the inverse Gaussian distribution
  */
  real inv_gaussian_rng(real mu, real lambda) {
    real z;
    real y;
    real x;
    real mu2;
    mu2 = square(mu);
    y = square(normal_rng(0,1));
    z = uniform_rng(0,1);
    x = mu + ( mu2 * y - mu * sqrt(4 * mu * lambda * y + mu2 * square(y)) )
      / (2 * lambda);
    if (z <= (mu / (mu + x))) return x;
    else return mu2 / x;
  }

  /** 
  * test function for csr_matrix_times_vector
  *
  * @param m Integer number of rows
  * @param n Integer number of columns
  * @param w Vector (see reference manual)
  * @param v Integer array (see reference manual)
  * @param u Integer array (see reference manual)
  * @param b Vector that is multiplied from the left by the CSR matrix
  * @return A vector that is the product of the CSR matrix and b
  */
  vector test_csr_matrix_times_vector(int m, int n, vector w, 
                                      int[] v, int[] u, vector b) {
    return csr_matrix_times_vector(m, n, w, v, u, b); 
  }

  /** 
   * Apply inverse link function to linear predictor
   * see help(binom) in R
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_bern(vector eta, int link) {
    vector[rows(eta)] pi;
    if (link < 1 || link > 5) 
      reject("Invalid link");
      
    if (link == 1)  // logit
      for(n in 1:rows(eta)) pi[n] = inv_logit(eta[n]);
    else if (link == 2)  // probit
      for(n in 1:rows(eta)) pi[n] = Phi(eta[n]);
    else if (link == 3)  // cauchit
      for(n in 1:rows(eta)) pi[n] = cauchy_cdf(eta[n], 0.0, 1.0);
    else if (link == 4)  // log
      for(n in 1:rows(eta)) pi[n] = exp(eta[n]);
    else if (link == 5)  // cloglog
      for(n in 1:rows(eta)) pi[n] = inv_cloglog(eta[n]);
    return pi;
  }

  /**
   * Increment with the unweighted log-likelihood
   * @param link An integer indicating the link function
   * @param eta0 A vector of linear predictors | y = 0
   * @param eta1 A vector of linear predictors | y = 1
   * @param N An integer array of length 2 giving the number of 
   *   observations where y = 0 and y = 1 respectively
   * @return lp__
   */
  real ll_bern_lp(vector eta0, vector eta1, int link, int[] N) {
    if (link < 1 || link > 5) 
      reject("Invalid link");
      
    if (link == 1) { // logit
      target += logistic_lccdf(eta0 | 0, 1);
      target += logistic_lcdf( eta1 | 0, 1);
    }
    else if (link == 2) {  // probit
      target += normal_lccdf(eta0 | 0, 1);
      target += normal_lcdf( eta1 | 0, 1);
    }
    else if (link == 3) {  // cauchit
      target += cauchy_lccdf(eta0 | 0, 1);
      target += cauchy_lcdf( eta1 | 0, 1);
    }
    else if(link == 4) {  // log
      vector[N[1]]       log_pi0;
      for (n in 1:N[1])  log_pi0[n] = log1m_exp(eta0[n]);
      target += log_pi0;
      target += eta1;  // already in log form
    }
    else if(link == 5) {  // cloglog
      vector[N[2]]       log_pi1;
      for (n in 1:N[2])  log_pi1[n] = log1m_exp(-exp(eta1[n]));
      target += log_pi1;
      target += -exp(eta0);
    }
    return target();
  }

  /** 
   * Pointwise (pw) log-likelihood vector
   *
   * @param y The integer outcome variable. Note that function is
   *  called separately with y = 0 and y = 1
   * @param eta Vector of linear predictions
   * @param link An integer indicating the link function
   * @return A vector
   */
  vector pw_bern(int y, vector eta, int link) {
    vector[rows(eta)] ll;
    if (link < 1 || link > 5) 
      reject("Invalid link");
      
    if (link == 1) {  // logit
      for (n in 1:rows(eta)) ll[n] = bernoulli_logit_lpmf(y | eta[n]);
    }
    else {  // link = probit, cauchit, log, or cloglog 
            // Note: this may not be numerically stable
      vector[rows(eta)] pi;
      pi = linkinv_bern(eta, link);
      for (n in 1:rows(eta)) ll[n] = bernoulli_lpmf(y | pi[n]);
    }
    return ll;
  }
  
  /** 
   * Apply inverse link function to linear predictor
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_binom(vector eta, int link) {
    vector[rows(eta)] pi;
    if (link < 1 || link > 5) 
      reject("Invalid link");
      
    if (link == 1)  // logit
      for(n in 1:rows(eta)) pi[n] = inv_logit(eta[n]);
    else if (link == 2)  // probit
      for(n in 1:rows(eta)) pi[n] = Phi(eta[n]);
    else if (link == 3)  // cauchit
      for(n in 1:rows(eta)) pi[n] = cauchy_cdf(eta[n], 0.0, 1.0);
    else if (link == 4)  // log 
      for(n in 1:rows(eta)) pi[n] = exp(eta[n]);
    else if (link == 5)  // cloglog
      for(n in 1:rows(eta)) pi[n] = inv_cloglog(eta[n]);
    return pi;
  }
  
  /**
  * Increment with the unweighted log-likelihood
  * @param y An integer array indicating the number of successes
  * @param trials An integer array indicating the number of trials
  * @param eta A vector of linear predictors
  * @param link An integer indicating the link function
  * @return lp__
  */
  real ll_binom_lp(int[] y, int[] trials, vector eta, int link) {
    if (link < 1 || link > 5) 
      reject("Invalid link");
      
    if (link == 1) target += binomial_logit_lpmf(y | trials, eta);
    else if (link <  4) target += binomial_lpmf( y | trials, linkinv_binom(eta, link));
    else if (link == 4) {  // log
      for (n in 1:num_elements(y)) {
        target += y[n] * eta[n];
        target += (trials[n] - y[n]) * log1m_exp(eta[n]);
      }
    }
    else if (link == 5) {  // cloglog
      real neg_exp_eta;
      for (n in 1:num_elements(y)) {
        neg_exp_eta = -exp(eta[n]);
        target += y[n] * log1m_exp(neg_exp_eta);
        target += (trials[n] - y[n]) * neg_exp_eta;
      }
    }
    return target();
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_binom(int[] y, int[] trials, vector eta, int link) {
    vector[rows(eta)] ll;
    if (link < 1 || link > 5) reject("Invalid link");
    if (link == 1) {  // logit
      for (n in 1:rows(eta)) 
        ll[n] = binomial_logit_lpmf(y[n] | trials[n], eta[n]);
    }
    else {  // link = probit, cauchit, log, or cloglog (unstable)
      vector[rows(eta)] pi;
      pi = linkinv_binom(eta, link);
      for (n in 1:rows(eta)) ll[n] = binomial_lpmf(y[n] | trials[n], pi[n]) ;
    }
    return ll;
  }

  vector linkinv_count(vector eta, int link) {
    vector[rows(eta)] phi;
    if (link < 1 || link > 3) 
      reject("Invalid link");
      
    if (link == 1) return exp(eta);  // log
    else if (link == 2) return eta;  // identity
    else  // link = sqrt
      for (n in 1:rows(eta)) phi[n] = square(eta[n]); 
    return phi;
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector for the Poisson distribution
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_pois(int[] y, vector eta, int link) {
    vector[rows(eta)] ll;
    if (link < 1 || link > 3) 
      reject("Invalid link");
      
    if (link == 1)  // log
      for (n in 1:rows(eta)) ll[n] = poisson_log_lpmf(y[n] | eta[n]);
    else {  // link = identity or sqrt
      vector[rows(eta)] phi;
      phi = linkinv_count(eta, link);
      for (n in 1:rows(eta)) ll[n] = poisson_lpmf(y[n] | phi[n]) ;
    }
    return ll;
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector for the negative binomial  distribution
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_nb(int[] y, vector eta, real theta, int link) {
    vector[rows(eta)] ll;
    vector[rows(eta)] rho;
    if (link < 1 || link > 3) 
      reject("Invalid link");
      
    rho = linkinv_count(eta, link);
    for (n in 1:rows(eta)) ll[n] = neg_binomial_2_lpmf(y[n] | rho[n], theta);
    return ll;
  }
  
  /** 
   * Create group-specific block-diagonal Cholesky factor, see section 2 of
   * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
   * @param len_theta_L An integer indicating the length of returned vector, 
   *   which lme4 denotes as m
   * @param p A matrix with the number variables on the LHS of each | in each
   *   longitudinal submodel. The rows correspond to each |, meaning the separate
   *   equations for each grouping variable, and the columns correspond to each
   *   longitudinal submodel. If subject ID is the only grouping variable then the
   *   matrix will have one row. If the joint model only has one longitudinal 
   *   submodel then the matrix will have one column.
   * @param dispersion Scalar standard deviation of the errors, called sigma by lme4
   * @param tau Vector of scale parameters whose squares are proportional to the 
   *   traces of the relative covariance matrices of the group-specific terms
   * @param scale Vector of prior scales that are multiplied by elements of tau
   * @param zeta Vector of positive parameters that are normalized into simplexes
   *   and multiplied by the trace of the covariance matrix to produce variances
   * @param rho Vector of radii in the onion method for creating Cholesky factors
   * @param z_T Vector used in the onion method for creating Cholesky factors
   * @return A vector that corresponds to theta in lme4
   */
  vector make_theta_L(int len_theta_L, int[] p, real dispersion,
                      vector tau, vector scale, vector zeta,
                      vector rho, vector z_T) {
    vector[len_theta_L] theta_L;
    int zeta_mark;
    int rho_mark;
    int z_T_mark;
    int theta_L_mark;
    zeta_mark = 1;
    rho_mark = 1;
    z_T_mark = 1;
    theta_L_mark = 1;
    
    // each of these is a diagonal block of the implicit Cholesky factor
    for (i in 1:size(p)) { 
      int nc;
      nc = p[i];
      if (nc == 1) { // "block" is just a standard deviation
        theta_L[theta_L_mark] = tau[i] * scale[i] * dispersion;
        // unlike lme4, theta[theta_L_mark] includes the dispersion term in it
        theta_L_mark = theta_L_mark + 1;
      }
      else { // block is lower-triangular               
        matrix[nc,nc] T_i; 
        real trace_T_i;
        vector[nc] pi; // variance = proportion of trace_T_i
        real std_dev;
        real T21;
        
        trace_T_i = square(tau[i] * scale[i] * dispersion) * nc;
        // unlike lme4, T_i includes the dispersion term in it
        pi = segment(zeta, zeta_mark, nc); // zeta ~ gamma(shape, 1)
        pi = pi / sum(pi);                 // thus pi ~ dirichlet(shape)
        zeta_mark = zeta_mark + nc;
        std_dev = sqrt(pi[1] * trace_T_i);
        T_i[1,1] = std_dev;
        
        // Put a correlation into T_i[2,1] and scale by std_dev
        std_dev = sqrt(pi[2] * trace_T_i);
        T21 = 2.0 * rho[rho_mark] - 1.0;
        rho_mark = rho_mark + 1;
        T_i[2,2] = std_dev * sqrt(1.0 - square(T21));
        T_i[2,1] = std_dev * T21;
        
        for (r in 2:(nc - 1)) { // scaled onion method to fill T_i
          int rp1;
          vector[r] T_row;
          real scale_factor;
          T_row = segment(z_T, z_T_mark, r);
          z_T_mark = z_T_mark + r;
          rp1 = r + 1;
          std_dev = sqrt(pi[rp1] * trace_T_i);
          scale_factor = sqrt(rho[rho_mark] / dot_self(T_row)) * std_dev;
          for(c in 1:r) T_i[rp1,c] = T_row[c] * scale_factor;
          T_i[rp1,rp1] = sqrt(1.0 - rho[rho_mark]) * std_dev;
          rho_mark = rho_mark + 1;
        }
        
        // now vech T_i
        for (c in 1:nc) for (r in c:nc) {
          theta_L[theta_L_mark] = T_i[r,c];
          theta_L_mark = theta_L_mark + 1;
        }
      }
    }
    return theta_L;
  }
 
  /** 
  * Create group-specific coefficients, see section 2 of
  * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  *
  * @param z_b Vector whose elements are iid normal(0,sigma) a priori
  * @param theta Vector with covariance parameters as defined in lme4
  * @param p An integer array with the number variables on the LHS of each |
  * @param l An integer array with the number of levels for the factor(s) on 
  *   the RHS of each |
  * @return A vector of group-specific coefficients
  */
  vector make_b(vector z_b, vector theta_L, int[] p, int[] l) {
    vector[rows(z_b)] b;
    int b_mark;
    int theta_L_mark;
    b_mark = 1;
    theta_L_mark = 1;
    for (i in 1:size(p)) {
      int nc;
      nc = p[i];
      if (nc == 1) {
        real theta_L_start;
        theta_L_start = theta_L[theta_L_mark];
        for (s in b_mark:(b_mark + l[i] - 1)) 
          b[s] = theta_L_start * z_b[s];
        b_mark = b_mark + l[i];
        theta_L_mark = theta_L_mark + 1;
      }
      else {
        matrix[nc,nc] T_i;
        T_i = rep_matrix(0, nc, nc);
        for (c in 1:nc) {
          T_i[c,c] = theta_L[theta_L_mark];
          theta_L_mark = theta_L_mark + 1;
          for(r in (c+1):nc) {
            T_i[r,c] = theta_L[theta_L_mark];
            theta_L_mark = theta_L_mark + 1;
          }
        }
        for (j in 1:l[i]) {
          vector[nc] temp;
          temp = T_i * segment(z_b, b_mark, nc);
          b_mark = b_mark - 1;
          for (s in 1:nc) b[b_mark + s] = temp[s];
          b_mark = b_mark + nc + 1;
        }
      }
    }
    return b;
  }

  /** 
   * Prior on group-specific parameters
   *
   * @param z_b A vector of primitive coefficients
   * @param z_T A vector of primitives for the unit vectors in the onion method
   * @param rho A vector radii for the onion method
   * @param zeta A vector of primitives for the simplexes
   * @param tau A vector of scale parameters
   * @param regularization A real array of LKJ hyperparameters
   * @param delta A real array of concentration paramters
   * @param shape A vector of shape parameters
   * @param t An integer indicating the number of group-specific terms
   * @param p An integer array with the number variables on the LHS of each |
   * @return nothing
   */
  void decov_lp(vector z_b, vector z_T, vector rho, vector zeta, vector tau,
                real[] regularization, real[] delta, vector shape,
                int t, int[] p) {
    int pos_reg;
    int pos_rho;
    target += normal_lpdf(z_b | 0, 1);
    target += normal_lpdf(z_T | 0, 1);
    pos_reg = 1;
    pos_rho = 1;
    for (i in 1:t) if (p[i] > 1) {
      vector[p[i] - 1] shape1;
      vector[p[i] - 1] shape2;
      real nu;
      nu = regularization[pos_reg] + 0.5 * (p[i] - 2);
      pos_reg = pos_reg + 1;
      shape1[1] = nu;
      shape2[1] = nu;
      for (j in 2:(p[i]-1)) {
        nu = nu - 0.5;
        shape1[j] = 0.5 * j;
        shape2[j] = nu;
      }
      target += beta_lpdf(rho[pos_rho:(pos_rho + p[i] - 2)] | shape1, shape2);
      pos_rho = pos_rho + p[i] - 1;
    }
    target += gamma_lpdf(zeta | delta, 1);
    target += gamma_lpdf(tau  | shape, 1);
  }
  
  /** 
   * Elementwise square root
   *
   * @param y A vector of non-negative numbers
   * @return A vector of square roots
   */
  vector sqrt_vec(vector y) {
    vector[rows(y)] out;
    for (i in 1:rows(y)) out[i] = sqrt(out[i]);
    return out;
  }

  /** 
   * Hierarchical shrinkage parameterization
   *
   * @param z_beta A vector of primitive coefficients
   * @param global A real array of positive numbers
   * @param local A vector array of positive numbers
   * @return A vector of coefficientes
   */
  vector hs_prior(vector z_beta, real[] global, vector[] local) {
    vector[rows(z_beta)] lambda;
    int K;
    K = rows(z_beta);
    for (k in 1:K) lambda[k] = local[1][k] * sqrt(local[2][k]);
    return z_beta .* lambda * global[1] * sqrt(global[2]);
  }

  /** 
   * Hierarchical shrinkage plus parameterization
   *
   * @param z_beta A vector of primitive coefficients
   * @param global A real array of positive numbers
   * @param local A vector array of positive numbers
   * @return A vector of coefficientes
   */
  vector hsplus_prior(vector z_beta, real[] global, vector[] local) {
    vector[rows(z_beta)] lambda;
    vector[rows(z_beta)] lambda_plus;
    int K;
    K = rows(z_beta);
    for (k in 1:K) {
      lambda[k] = local[1][k] * sqrt(local[2][k]);
      lambda_plus[k] = local[3][k] * sqrt(local[4][k]);
    }
    return z_beta .* lambda .* lambda_plus * global[1] * sqrt(global[2]);
  }
  
  /** 
   * Divide a scalar by a vector
   *
   * @param x The scalar in the numerator
   * @param y The vector in the denominator
   * @return An elementwise vector
   */
  vector divide_real_by_vector(real x, vector y) {
    vector[rows(y)] ret;
    for (n in 1:rows(y)) ret[n] = x / y[n];
    return ret;
  }
 
  /** 
  * Reorder the vector of group-specific coefficients
  *
  * @param b Vector whose elements are ordered in the following nested way:
  *   factor level (e.g. "subject ID" or "clinic ID") within grouping 
  *   factor (e.g. "subjects" or "clinics")
  * @param p An integer array with the number variables on the LHS of each |
  * @param pmat A matrix with the number variables on the LHS of each | in each
  *   longitudinal submodel. The rows correspond to each |, meaning the separate
  *   equations for each grouping variable, and the columns correspond to each
  *   longitudinal submodel. If subject ID is the only grouping variable then the
  *   matrix will have one row. If the joint model only has one longitudinal 
  *   submodel then the matrix will have one column.
  * @param p An integer array with the number of random coefficients for the 
  *   LHS of each |
  * @param qmat A matrix with the number of random coefficients for the LHS of each | 
  *   The rows correspond to each |, meaning the separate equations for 
  *   each grouping variable, and the columns correspond to each longitudinal
  *   submodel.
  * @param l An integer array with the number of levels for the factor(s) on 
  *   the RHS of each |
  * @param M An integer specifying the number of longitudinal submodels
  * @return A vector (containing the group-specific coefficients) whose whose 
  *   elements are ordered in the following nested way: factor level, within
  *   grouping factor, within longitudinal submodel
  */ 
  vector reorder_b(vector b, int[] p, int[,] pmat, int[] q, 
                   int[,] qmat, int[] l, int M) {
    vector[rows(b)] b_new;
    int nq_mark;
    nq_mark = 1;
    for (i in 1:size(p)) { // loop over grouping factors 
      int np;  // num. of random effects for grouping factor
      int nq;  // num. of random coefs for grouping factor
      vector[q[i]] b_i;
      np = p[i];
      nq = q[i];
      b_i = segment(b, nq_mark, nq);
      nq_mark = nq_mark + nq;
      for (m in 1:M) { // loop over markers
        int shift_collect;
        int shift_store;
        int start_collect;
        int start_store;
        int end_store;
        if (m == 1) shift_collect = 0;
        else shift_collect = sum(pmat[i, 1:(m-1)]);
        if (m == 1) shift_store = 0;
        else shift_store = sum(qmat[i, 1:(m-1)]);    
        for (j in 1:l[i]) {
        start_collect = (j - 1) * pmat[i,m] + shift_collect + 1;
        start_store   = (j - 1) * pmat[i,m] + shift_store + 1;
        end_store     = (j - 1) * pmat[i,m] + shift_store + pmat[i,m];        
        b_new[start_store:end_store] = segment(b_i, start_collect, pmat[i,m]);
        }
      }
    }  
    return b_new;
  }
  
  /** 
  * Create a design matrix for a shared random effects association
  * structure in the joint model
  *
  * @param b Vector of group-specific coefficients
  * @param l An integer array with the number of levels for the factor(s) on 
  *   the RHS of each |
  * @param p An integer array with the number of variables on the LHS of each |
  * @param pmat A matrix with the number variables on the LHS of each | in each
  *   longitudinal submodel. The rows correspond to each |, meaning the separate
  *   equations for each grouping variable, and the columns correspond to each
  *   longitudinal submodel. If subject ID is the only grouping variable then the
  *   matrix will have one row. If the joint model only has one longitudinal 
  *   submodel then the matrix will have one column.
  * @param Npat Integer specifying number of individuals represented 
  *   in vector b
  * @param quadnodes The number of quadrature nodes
  * @param which_b Integer array specifying the indices
  *   of the random effects to use in the association structure
  * @param sum_size_which_b Integer specifying total number of 
  *   random effects that are to be used in the association structure
  * @param size_which_b Integer array specifying number of random effects from
  *   each long submodel that are to be used in the association structure
  * @param t_i Integer specifying the index of the grouping factor that
  *   corresponds to the patient-level
  * @param M An integer specifying the number of longitudinal submodels
  * @return A matrix with the desired random effects represented
  *   in columns, and the individuals on the rows; the matrix is
  *   repeated (quadnodes + 1) times (bounded by rows)
  */  
  matrix make_x_assoc_sh(vector b, int[] l, int[] p, int[,] pmat, int Npat, int quadnodes,
	                  int[] which_b, int sum_size_which_b, int[] size_which_b, 
                    int t_i, int M) {
    int prior_shift;    // num. ranefs prior to subject-specific ranefs
    int start_store;
    int end_store;	
    matrix[Npat,sum_size_which_b] temp;
    matrix[(Npat*(quadnodes+1)),sum_size_which_b] x_assoc_sh;						  
    if (t_i == 1) prior_shift = 0;
    else prior_shift = sum(l[1:(t_i-1)]);
    for (i in 1:Npat) {
      int mark;
      int start_collect;  // index start of subject-specific ranefs for patient
      mark = 1;
      start_collect = prior_shift + (i - 1) * p[t_i];
      for (m in 1:M) {
        if (size_which_b[m] > 0) {
          int shift;  // num. subject-specific ranefs in prior submodels
          int j_shift; // shift in indexing of which_b vector
          if (m == 1) {
            shift = 0;
            j_shift = 0;
          }
          else {
            shift = sum(pmat[t_i, 1:(m-1)]);
            j_shift = sum(size_which_b[1:(m-1)]);
          }
          for (j in (j_shift+1):(j_shift+size_which_b[m])) {
            int item_collect;   // subject-specific ranefs to select for current submodel
            item_collect = start_collect + shift + which_b[j];
            temp[i,mark] = b[item_collect];
            mark = mark + 1;
          }
        }      
      }
    }
    for (i in 1:(quadnodes+1)) {
      start_store = (i - 1) * Npat + 1;
      end_store   = i * Npat;		
      x_assoc_sh[start_store:end_store,] = temp;
    }
  return x_assoc_sh;
  }
}
data {
  // dimensions
  int<lower=1> M;  // num. of long. submodels
  int<lower=0> Npat;  // num. individuals (equal to l[1] - 1)
  int<lower=0> y_N[M];  // num. of obs. in each long. submodel
  int<lower=0> sum_y_N; // total num. of obs. across all long submodels
  int<lower=0> sum_y_real_N; // total num. of obs. across all long submodels with real outcomes
  int<lower=0> sum_y_int_N; // total num. of obs. across all long submodels with integer outcomes
  int<lower=0,upper=sum_y_N> y_beg[M]; // index of first obs. for each submodel
  int<lower=0,upper=sum_y_N> y_end[M]; // index of last obs. for each submodel
  int<lower=0,upper=sum_y_real_N> y_real_beg[M]; // index of first obs. for each submodel
  int<lower=0,upper=sum_y_real_N> y_real_end[M]; // index of last obs. for each submodel
  int<lower=0,upper=sum_y_int_N> y_int_beg[M]; // index of first obs. for each submodel
  int<lower=0,upper=sum_y_int_N> y_int_end[M]; // index of last obs. for each submodel
  int<lower=0> y_N01[M,2];  // num. of bernoulli 0/1 observations in each long. submodel
  int<lower=0> y_K[M];  // num. of predictors in each long. submodel
  int<lower=0> sum_y_K; // total num. of predictors across all long submodels
  int<lower=0> e_K;   // num. of predictors in event submodel
  int<lower=0> a_K;   // num. of association parameters
  int<lower=0> quadnodes;  // num. of nodes for Gauss-Kronrod quadrature 
  int<lower=0> Npat_times_quadnodes;
  int<lower=0,upper=M> sum_y_has_intercept; // num. submodels w/ intercept
  int<lower=0,upper=M> sum_y_has_intercept_unbound; // num. submodels w/ unbounded intercept
  int<lower=0,upper=M> sum_y_has_intercept_lobound; // num. submodels w/ lower bounded intercept
  int<lower=0,upper=M> sum_y_has_intercept_upbound; // num. submodels w/ upper bounded intercept
  int<lower=0,upper=M> sum_y_has_dispersion; // num. submodels w/ dispersion term
  
  // data for longitudinal submodel(s)
  int<lower=1> family[M];          // family
  int<lower=0,upper=1> any_fam_3;  // any long. submodel with family == 3
  int<lower=1> link[M];            // link function, varies by .stan file
  int<lower=0,upper=1> y_centre;         // 1 = yes for centred predictor matrix
  int<lower=0,upper=1> y_has_intercept[M];  // 1 = yes
  int<lower=0,upper=1> y_has_intercept_unbound[M];  // intercept unbounded
  int<lower=0,upper=1> y_has_intercept_lobound[M];  // intercept lower bound at 0
  int<lower=0,upper=1> y_has_intercept_upbound[M];  // intercept upper bound at 0
  int<lower=0,upper=1> y_has_weights[M];    // 1 = Yes
  //int<lower=0,upper=1> y_has_offset;     // 1 = Yes
  int<lower=0,upper=1> y_has_dispersion[M];    // 1 = Yes
  vector[sum_y_real_N] y_real;           // outcome vector, reals                
  int<lower=0> y_int[sum_y_int_N];        // outcome vector, integers                
  vector[sum_y_K*(y_centre>0)] y_xbar;  // predictor means
  matrix[sum_y_N,sum_y_K] y_X;   // predictor matrix, possibly centred          
  vector[sum_y_N] y_weights;  // weights, set to zero if not used
  int<lower=0> trials[sum_y_N];  // num. binomial trials, set to zero if not used
  //vector[sum_y_N*y_has_offset] y_offset; 
  int<lower=0> num_non_zero;   // number of non-zero elements in the Z matrix
  vector[num_non_zero] w;  // non-zero elements in the implicit Z matrix
  int<lower=0> v[num_non_zero]; // column indices for w  
  int<lower=0> u[(sum_y_N+1)];  // where the non-zeros start in each row 
  
  // data for event submodel
  int<lower=0,upper=1> basehaz_weibull;  // weibull baseline hazard
  int<lower=0,upper=1> basehaz_piecewise;  // piecewise constant baseline hazard
  int<lower=0,upper=1> basehaz_splines;  // cubic splines baseline hazard
  int<lower=0> splines_df;  // df for cubic splines baseline hazard
  int<lower=0,upper=1> e_centre;  // 1 = yes for centred predictor matrix
  int<lower=0,upper=1> e_has_intercept;  // 1 = yes
  int<lower=0> nrow_y_Xq;     // num. rows in long. predictor matrix at quad points
  int<lower=0> nrow_e_Xq;   // num. rows in event predictor matrix at quad points
  matrix[(M*nrow_y_Xq),sum_y_K] y_Xq; // predictor matrix (long submodel) at quadpoints, possibly centred              
  matrix[nrow_e_Xq,e_K] e_Xq;         // predictor matrix (event submodel) at quadpoints, possibly centred
  vector[nrow_e_Xq] e_times;          // event times and unstandardised quadrature points
  matrix[(nrow_e_Xq*basehaz_splines),splines_df] e_ns_times; // basis for cubic splines baseline hazard
  vector[nrow_e_Xq] e_d;              // event indicator, followed by dummy indicator for quadpoints
  vector[e_K*(e_centre>0)] e_xbar;   // predictor means (event submodel)
  int<lower=0> num_non_zero_Zq;    // number of non-zero elements in the Z matrix (at quadpoints)
  vector[num_non_zero_Zq] w_Zq;  // non-zero elements in the implicit Z matrix (at quadpoints)
  int<lower=0> v_Zq[num_non_zero_Zq]; // column indices for w (at quadpoints)
  int<lower=0> u_Zq[(M*nrow_y_Xq+1)]; // where the non-zeros start in each row (at quadpoints)
  vector[Npat_times_quadnodes] quadweight_times_half_eventtime;
    
  // data for association structure
  int<lower=0,upper=1> assoc;              // 0 = no jm association structure, 1 = any jm association structure
  int<lower=0,upper=1> has_assoc_ev[M];    // eta value
  int<lower=0,upper=1> has_assoc_es[M];    // eta slope
  int<lower=0,upper=1> has_assoc_cv[M];    // mu value
  int<lower=0,upper=1> has_assoc_cs[M];    // mu slope
  int<lower=0,upper=M> sum_has_assoc_ev;   // num. long submodels linked via eta value
  int<lower=0,upper=M> sum_has_assoc_es;   // num. long submodels linked via eta slope
  int<lower=0,upper=M> sum_has_assoc_cv;   // num. long submodels linked via mu value
  int<lower=0,upper=M> sum_has_assoc_cs;   // num. long submodels linked via mu slope 
  int<lower=0> sum_size_which_b;           // num. of shared random effects
  int<lower=0> size_which_b[M];            // num. of shared random effects for each long submodel
  int<lower=1> which_b[sum_size_which_b];  // which random effects are shared for each long submodel

  // data for calculating slope
  real<lower=0> eps;  // time shift used for numerically calculating derivative
  matrix[(M*nrow_y_Xq*((sum_has_assoc_es + sum_has_assoc_cs) > 0)),sum_y_K] 
    y_Xq_eps; // predictor matrix (long submodel) at quadpoints plus time shift of epsilon              
  int<lower=0> num_non_zero_Zq_eps;        // number of non-zero elements in the Zq_eps matrix (at quadpoints plus time shift of epsilon)
  vector[num_non_zero_Zq_eps] w_Zq_eps;    // non-zero elements in the implicit Zq_eps matrix (at quadpoints plus time shift of epsilon)
  int<lower=0> v_Zq_eps[num_non_zero_Zq_eps]; // column indices for w (at quadpoints plus time shift of epsilon)
  int<lower=0> u_Zq_eps[(M*nrow_y_Xq*((sum_has_assoc_es + sum_has_assoc_cs) > 0) + 1)]; 
    // where the non-zeros start in each row (at quadpoints plus time shift of epsilon)

  // data for random effects model
  int<lower=1> t;     	        // num. of grouping factors
  int<lower=1,upper=t> t_i;     // index of grouping factor corresponding to patient-level
  int<lower=0> pmat[t,M];       // num. random effects for each grouping factor (t) in each submodel (M)
  int<lower=0> p[t];            // total num. random effects for each grouping factor (t) (rowsums of pmat)
  int<lower=1> l[t];            // num. levels for each grouping factor
  int<lower=0> qmat[t,M];       // = l * pmat --> num. random coefs for each grouping factor in each submodel
  int<lower=0> q[t];            // = l * p --> num. random coefs for each grouping factor
  int<lower=0> len_theta_L;     // length of the theta_L vector
  int<lower=0> len_b;           // length of the b vector

  // priors: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus
  int<lower=0,upper=4> priorLong_dist;
  int<lower=0,upper=2> priorLong_dist_for_intercept;  
  int<lower=0,upper=4> priorEvent_dist;
  int<lower=0,upper=2> priorEvent_dist_for_intercept;
  int<lower=0,upper=4> priorAssoc_dist;
    
  // hyperparameters for priors, set to 0 if there is no prior
  vector[sum_y_K] priorLong_mean;
  vector[M]     priorLong_mean_for_intercept;
  vector[e_K]   priorEvent_mean;
  real          priorEvent_mean_for_intercept;
  vector[a_K]   priorAssoc_mean;
  vector<lower=0>[sum_y_K] priorLong_scale;
  vector<lower=0>[M]     priorLong_scale_for_intercept;
  vector<lower=0>[e_K]   priorEvent_scale;
  real<lower=0>          priorEvent_scale_for_intercept;
  vector<lower=0>[a_K]   priorAssoc_scale;
  vector<lower=0>[sum_y_K] priorLong_df;
  vector<lower=0>[M]     priorLong_df_for_intercept;
  vector<lower=0>[e_K]   priorEvent_df;
  real<lower=0>          priorEvent_df_for_intercept; 
  vector<lower=0>[a_K]   priorAssoc_df;
  vector<lower=0>[sum_y_has_dispersion] priorLong_scale_for_dispersion;
  real<lower=0> priorEvent_scale_for_weibull[basehaz_weibull];
  vector<lower=0>[splines_df] priorEvent_scale_for_splines;
 
  // hyperparameters for random effects model
  vector<lower=0>[t] shape; 
  vector<lower=0>[t] scale;
  int<lower=0> len_concentration;
  real<lower=0> concentration[len_concentration];
  int<lower=0> len_regularization;
  real<lower=0> regularization[len_regularization];
   
  // flag indicating whether to draw from the prior
  int<lower=0,upper=1> prior_PD;  // 1 = yes

}
transformed data {
  real poisson_max;
  vector[sum_y_real_N] sqrt_y;
  vector[sum_y_real_N] log_y;
  real sum_log_y[M];
  vector[nrow_e_Xq] e_log_times;  // log of event times and unstandardised quadrature points
  int<lower=0,upper=1> y_t_any_124;
  int<lower=0,upper=1> y_t_all_124;
  int<lower=0,upper=1> e_t_any_124;  
  int<lower=0,upper=1> e_t_all_124;   
  int<lower=0,upper=1> a_t_any_124;  
  int<lower=0,upper=1> a_t_all_124; 
  int<lower=0> y_hs;
  int<lower=0> e_hs;                 
  int<lower=0> a_hs;                 
  int<lower=0> len_z_T;
  int<lower=0> len_var_group;
  int<lower=0> len_rho;
  real<lower=0> delta[len_concentration];
  int<lower=1> pos;
 
  poisson_max = pow(2.0, 30.0); 
  
  // calculate transformations of outcome
  for (m in 1:M) {
    if (family[m] == 1) {
      for (n in y_real_beg[m]:y_real_end[m]) {
		    sqrt_y[n] = not_a_number();
		    log_y[n] = not_a_number();
	    }
      sum_log_y[m] = not_a_number();
	  }
    else if (family[m] == 2) {
      for (n in y_real_beg[m]:y_real_end[m]) {
		    sqrt_y[n] = not_a_number();
		    log_y[n] = log(y_real[n]);
	    }
      sum_log_y[m] = sum(log_y);
    }
    else if (family[m] == 3) {
      for (n in y_real_beg[m]:y_real_end[m]) {
		    sqrt_y[n] = sqrt(y_real[n]);
		    log_y[n] = log(y_real[n]);
	     }
      sum_log_y[m] = sum(log_y);
    }
    else sum_log_y[m] = not_a_number();
  }
  
  // calculate log of event times and unstandardised quadpoints
  e_log_times = log(e_times);   
 
  // priors for longitudinal submodels  
  if (priorLong_dist <= 2) y_hs = 0;
  else if (priorLong_dist == 3) y_hs = 2;
  else if (priorLong_dist == 4) y_hs = 4;
  if (priorLong_dist == 2) {
    y_t_any_124 = 0;
    y_t_all_124 = 1;
    for (k in 1:sum_y_K) {
      if (priorLong_df[k] == 1 || priorLong_df[k] == 2 || priorLong_df[k] == 4)
        y_t_any_124 = 1;
      else y_t_all_124 = 0;
    }
  }
  else {
    y_t_any_124 = 0;
    y_t_all_124 = 0;
  }

  // priors for event submodel
  if (priorEvent_dist <= 2) e_hs = 0;
  else if (priorEvent_dist == 3) e_hs = 2;
  else if (priorEvent_dist == 4) e_hs = 4;   
  if (priorEvent_dist == 2) {
    e_t_any_124 = 0;
    e_t_all_124 = 1;
    for (k in 1:e_K) {
      if (priorEvent_df[k] == 1 || priorEvent_df[k] == 2 || priorEvent_df[k] == 4)
        e_t_any_124 = 1;
      else e_t_all_124 = 0;
    }
  }
  else {
    e_t_any_124 = 0;
    e_t_all_124 = 0;
  }
  
  // priors for association parameters
  if (priorAssoc_dist <= 2) a_hs = 0;
  else if (priorAssoc_dist == 3) a_hs = 2;
  else if (priorAssoc_dist == 4) a_hs = 4;   
  if (priorAssoc_dist == 2) {
    a_t_any_124 = 0;
    a_t_all_124 = 1;
    for (k in 1:a_K) {
      if (priorAssoc_df[k] == 1 || priorAssoc_df[k] == 2 || priorAssoc_df[k] == 4)
        a_t_any_124 = 1;
      else a_t_all_124 = 0;
    }
  }
  else {
    a_t_any_124 = 0;
    a_t_all_124 = 0;
  }  
  
  // prior for covariance
  len_z_T = 0;
  len_var_group = sum(p) * (t > 0);
  len_rho = sum(p) - t;
  pos = 1;
  for (i in 1:t) {
    if (p[i] > 1) {
      for (j in 1:p[i]) {
        delta[pos] = concentration[j];
        pos = pos + 1;
      }
    }
    for (j in 3:p[i]) len_z_T = len_z_T + p[i] - 1;
  }
  
}
parameters {
  // parameters for longitudinal submodel(s)
  real y_gamma_unbound[sum_y_has_intercept_unbound];
  real<lower=0> y_gamma_lobound[sum_y_has_intercept_lobound];  
  real<upper=0> y_gamma_upbound[sum_y_has_intercept_upbound];  

  vector[sum_y_K] y_z_beta;                // primative coefs (long submodels)
  vector<lower=0>[sum_y_has_dispersion] y_dispersion_unscaled; # interpretation depends on family!
  #vector<lower=0>[sum_y_noise_N] y_noise; // do not store this
  
  // parameters for event submodel
  real e_gamma[e_has_intercept];          // intercept (event model)
  vector[e_K] e_z_beta;                   // primative coefs (event submodel)
  real<lower=0> weibull_shape_unscaled[basehaz_weibull];  // unscaled weibull shape parameter 
  vector[splines_df] splines_coefs_unscaled;       // unscaled coefs for cubic splines baseline hazard
 
  // parameters for association structure
  vector[a_K] a_z_beta;   // primative coefs
    
  // parameters for random effects model
  vector[len_b] z_b;
  vector[len_z_T] z_T;
  vector<lower=0,upper=1>[len_rho] rho;
  vector<lower=0>[len_concentration] zeta;
  vector<lower=0>[t] tau;  

  // parameters for priors
  real<lower=0> y_global[y_hs];
  vector<lower=0>[(y_hs>0)*sum_y_K] y_local[y_hs];
  real<lower=0> e_global[e_hs];
  vector<lower=0>[(e_hs>0)*e_K] e_local[e_hs];
  real<lower=0> a_global[a_hs];
  vector<lower=0>[(a_hs>0)*a_K] a_local[a_hs];
  
}
transformed parameters {
  // parameters for longitudinal submodel(s)
  vector[sum_y_K] y_beta;
  vector[sum_y_has_dispersion] y_dispersion;
  
  // parameters for event submodel
  vector[e_K] e_beta; 
  real weibull_shape[basehaz_weibull];
  vector[splines_df] splines_coefs;     
  
  // parameters for GK quadrature  
  vector[(M*nrow_y_Xq)] y_eta_q;          // linear predictor (all long submodels) evaluated at quadpoints
  vector[(M*nrow_y_Xq)*((sum_has_assoc_es + sum_has_assoc_cs) > 0)] y_eta_q_eps; 
    // linear predictor (all long submodels) evaluated at quadpoints plus time shift of epsilon
  vector[nrow_e_Xq] e_eta_q;      // linear predictor (event submodel) evaluated at quadpoints
  vector[nrow_e_Xq] log_basehaz;      // baseline hazard evaluated at quadpoints
  vector[nrow_e_Xq] ll_haz_q;     // log hazard contribution to the log likelihood for the event model at event time and quad points
  vector[Npat] ll_haz_eventtime;  // log hazard contribution to the log likelihood for the event model AT the event time only
  vector[Npat_times_quadnodes] ll_haz_quadtime;    // log hazard for the event model AT the quadrature points only
  vector[Npat_times_quadnodes] ll_surv_eventtime;  // log survival contribution to the log likelihood for the event model AT the event time
  real sum_ll_haz_eventtime;
  real sum_ll_surv_eventtime;
  real ll_event;                                   // log likelihood for the event model    
  
  // parameters for association structure  
  vector[a_K] a_beta;           

  // parameters for random effects model
  vector[len_theta_L] theta_L; 
  vector[len_b] b; 
  vector[len_b] b_by_model; 

  // parameters for longitudinal submodel(s)
  if      (priorLong_dist == 0) y_beta = y_z_beta;
  else if (priorLong_dist == 1) y_beta = y_z_beta .* priorLong_scale + priorLong_mean;
  else if (priorLong_dist == 2) for (k in 1:sum_y_K) {
    real P;
    if (priorLong_df[k] == 1) {
      P = Phi(y_z_beta[k]);
      y_beta[k] = tan(pi() * (P - 0.5));
    }
    else if (priorLong_df[k] == 2) {
      P = Phi(y_z_beta[k]);
      y_beta[k] = 2 * (P - 0.5) / sqrt(2.0 * P * (1 - P));
    }
    else if (priorLong_df[k] == 4) {
      real q_a;
      P = Phi(y_z_beta[k]);
      q_a = sqrt(4.0 * P * (1 - P));
      q_a = cos(acos(q_a) / 3) / q_a;
      y_beta[k] = 2 * sqrt(q_a - 1);
      if (P < 0.5) y_beta[k] = -y_beta[k];
    }
    else y_beta[k] = y_z_beta[k];
    y_beta[k] = y_beta[k] * priorLong_scale[k] + priorLong_mean[k];
  }
  else if (priorLong_dist == 3) y_beta = hs_prior(y_z_beta, y_global, y_local);
  else if (priorLong_dist == 4) y_beta = hsplus_prior(y_z_beta, y_global, y_local);
  
  if (sum_y_has_dispersion > 0) {
    int mark;
    mark = 1;
    for (m in 1:M) {
      if (y_has_dispersion[m] == 1) {
        if (priorLong_scale_for_dispersion[mark] > 0)
  	      y_dispersion[mark] =  priorLong_scale_for_dispersion[mark] * y_dispersion_unscaled[mark];
        else 
          y_dispersion[mark] = y_dispersion_unscaled[mark];
        mark = mark + 1;
      }    
    }    
  }
  
  // parameters for event submodel
  if      (priorEvent_dist == 0) e_beta = e_z_beta;
  else if (priorEvent_dist == 1) e_beta = e_z_beta .* priorEvent_scale + priorEvent_mean;
  else if (priorEvent_dist == 2) for (k in 1:e_K) {
    real P;
    if (priorEvent_df[k] == 1) {
      P = Phi(e_z_beta[k]);
      e_beta[k] = tan(pi() * (P - 0.5));
    }
    else if (priorEvent_df[k] == 2) {
      P = Phi(e_z_beta[k]);
      e_beta[k] = 2 * (P - 0.5) / sqrt(2.0 * P * (1 - P));
    }
    else if (priorEvent_df[k] == 4) {
      real q_a;
      P = Phi(e_z_beta[k]);
      q_a = sqrt(4.0 * P * (1 - P));
      q_a = cos(acos(q_a) / 3) / q_a;
      e_beta[k] = 2 * sqrt(q_a - 1);
      if (P < 0.5) e_beta[k] = -e_beta[k];
    }
    else e_beta[k] = e_z_beta[k];
    e_beta[k] = e_beta[k] * priorEvent_scale[k] + priorEvent_mean[k];
  }
  else if (priorEvent_dist == 3) e_beta = hs_prior(e_z_beta, e_global, e_local);
  else if (priorEvent_dist == 4) e_beta = hsplus_prior(e_z_beta, e_global, e_local);
  
  if (basehaz_weibull == 1) {
    if (priorEvent_scale_for_weibull[1] > 0)
      weibull_shape[1] = priorEvent_scale_for_weibull[1] * weibull_shape_unscaled[1];
    else weibull_shape[1] = weibull_shape_unscaled[1];
  } else if (basehaz_splines == 1) {
    splines_coefs = priorEvent_scale_for_splines .* splines_coefs_unscaled;
  }   
  
  // parameters for association structure
  if      (priorAssoc_dist == 0) a_beta = a_z_beta;
  else if (priorAssoc_dist == 1) a_beta = a_z_beta .* priorAssoc_scale + priorAssoc_mean;
  else if (priorAssoc_dist == 2) for (k in 1:a_K) {
    real P;
    if (priorAssoc_df[k] == 1) {
      P = Phi(a_z_beta[k]);
      a_beta[k] = tan(pi() * (P - 0.5));
    }
    else if (priorAssoc_df[k] == 2) {
      P = Phi(a_z_beta[k]);
      a_beta[k] = 2 * (P - 0.5) / sqrt(2.0 * P * (1 - P));
    }
    else if (priorAssoc_df[k] == 4) {
      real q_a;
      P = Phi(a_z_beta[k]);
      q_a = sqrt(4.0 * P * (1 - P));
      q_a = cos(acos(q_a) / 3) / q_a;
      a_beta[k] = 2 * sqrt(q_a - 1);
      if (P < 0.5) a_beta[k] = -a_beta[k];
    }
    else a_beta[k] = a_z_beta[k];
    a_beta[k] = a_beta[k] * priorAssoc_scale[k] + priorAssoc_mean[k];
  }
  else if (priorAssoc_dist == 3) a_beta = hs_prior(a_z_beta, a_global, a_local);
  else if (priorAssoc_dist == 4) a_beta = hsplus_prior(a_z_beta, a_global, a_local);

  // parameters for random effects model
  if (t > 0) {
    theta_L = make_theta_L(len_theta_L, p, 1.0, tau, scale, zeta, rho, z_T);
    b = make_b(z_b, theta_L, p, l);
    if (M > 1) b_by_model = reorder_b(b, p, pmat, q, qmat, l, M);
	  else b_by_model = b;
  }
  
  //===============
  // GK quadrature
  //===============
 
 
  // Longitudinal submodel(s): linear predictor at event and quad times
  if (sum_y_K > 0) y_eta_q = y_Xq * y_beta;
  else y_eta_q = rep_vector(0.0, (M*nrow_y_Xq));
  //if (y_has_offset == 1) y_eta_q = y_eta_q + y_offset;
  y_eta_q = y_eta_q + csr_matrix_times_vector((M*nrow_y_Xq), len_b, w_Zq, v_Zq, u_Zq, b_by_model);

  // Longitudinal submodel(s): slope of linear predictor at event and quad times
  if ((sum_has_assoc_es > 0) || (sum_has_assoc_cs > 0)) {
    if (sum_y_K > 0) y_eta_q_eps = y_Xq_eps * y_beta;
    else y_eta_q_eps = rep_vector(0.0, (M*nrow_y_Xq));
    // !!! if (y_has_offset == 1) y_eta_q_eps = y_eta_q_eps + y_offset; # ignore offset in derivative?
    y_eta_q_eps = y_eta_q_eps + csr_matrix_times_vector((M*nrow_y_Xq), len_b, w_Zq_eps, v_Zq_eps, u_Zq_eps, b_by_model);
  }
  // Event submodel: linear predictor at event and quad times
  if (e_K > 0) e_eta_q = e_Xq * e_beta;
  else e_eta_q = rep_vector(0.0, nrow_e_Xq);
  if (e_has_intercept == 1) {
    e_eta_q = e_eta_q + e_gamma[1];
  }
  else if (e_centre == 1) {
    // correction to eta if model has no intercept (because X is centered)
    e_eta_q = e_eta_q + dot_product(e_xbar, e_beta); 
  }
  
  if (assoc == 1) {
    int mark;
	  mark = 0;
    for (m in 1:M) {
      vector[nrow_y_Xq] ysep_q;         // expected long. outcome at event and quad times   
      vector[nrow_y_Xq] ysep_q_eps;     // expected long. outcome at event and quad times plus time shift of epsilon 
      
      # Prep work for association structures
      
      # Linear predictor
      if (y_has_intercept[m] == 1) {
        if (y_has_intercept_unbound[m] == 1) 
          y_eta_q = y_eta_q +
                             y_gamma_unbound[sum(y_has_intercept_unbound[1:m])];
        else if (y_has_intercept_lobound[m] == 1)
          y_eta_q = y_eta_q - min(y_eta_q) + 
                             y_gamma_lobound[sum(y_has_intercept_lobound[1:m])];
  	  else if (y_has_intercept_upbound[m] == 1)
          y_eta_q = y_eta_q - max(y_eta_q) + 
                             y_gamma_lobound[sum(y_has_intercept_lobound[1:m])];
      }
      else if (y_centre == 1) {
        int mark_beg;
        int mark_end;
        if (m == 1) mark_beg = 1;
        else mark_beg = sum(y_K[1:(m-1)]) + 1;
        mark_end = sum(y_K[1:m]);   
        // correction to eta if model has no intercept (and X is centered)
        y_eta_q = y_eta_q + 
              dot_product(y_xbar[mark_beg:mark_end], y_beta[mark_beg:mark_end]); 
      }
      
      # Linear predictor at time plus epsilon
      if ((has_assoc_es[m] == 1) || (has_assoc_cs[m] == 1)) {
        if (y_has_intercept[m] == 1) {
          if (y_has_intercept_unbound[m] == 1) 
            y_eta_q_eps = y_eta_q_eps +
                               y_gamma_unbound[sum(y_has_intercept_unbound[1:m])];
          else if (y_has_intercept_lobound[m] == 1)
            y_eta_q_eps = y_eta_q_eps - min(y_eta_q_eps) + 
                               y_gamma_lobound[sum(y_has_intercept_lobound[1:m])];
    	    else if (y_has_intercept_upbound[m] == 1)
            y_eta_q_eps = y_eta_q_eps - max(y_eta_q_eps) + 
                               y_gamma_lobound[sum(y_has_intercept_lobound[1:m])];
        }
        else if (y_centre == 1) {
          int mark_beg;
          int mark_end;
          if (m == 1) mark_beg = 1;
          else mark_beg = sum(y_K[1:(m-1)]) + 1;
          mark_end = sum(y_K[1:m]);   
          // correction to eta if model has no intercept (and X is centered)
          y_eta_q_eps = y_eta_q_eps + 
                dot_product(y_xbar[mark_beg:mark_end], y_beta[mark_beg:mark_end]); 
        }
      } 
      
      # Expected value 
      if ((has_assoc_cv[m] == 1) || (has_assoc_cs[m] == 1)) {
        if (family[m] == 1) 
          ysep_q = linkinv_gauss(y_eta_q, link[m]);
        else if (family[m] == 2) 
          ysep_q = linkinv_gamma(y_eta_q, link[m]);
        else if (family[m] == 3)
          ysep_q = linkinv_inv_gaussian(y_eta_q, link[m]);
        else if (family[m] == 4)
          ysep_q = linkinv_bern(y_eta_q, link[m]);	
        else if (family[m] == 5)		  
		      ysep_q = linkinv_binom(y_eta_q, link[m]); 
      }	      
      
      # Expected value at time plus epsilon
       if (has_assoc_cs[m] == 1) {
        if (family[m] == 1) 
          ysep_q_eps = linkinv_gauss(y_eta_q_eps, link[m]);
        else if (family[m] == 2) 
          ysep_q_eps = linkinv_gamma(y_eta_q_eps, link[m]);
        else if (family[m] == 3)
          ysep_q_eps = linkinv_inv_gaussian(y_eta_q_eps, link[m]);
        else if (family[m] == 4)
          ysep_q_eps = linkinv_bern(y_eta_q_eps, link[m]);	
        else if (family[m] == 5)		  
		      ysep_q_eps = linkinv_binom(y_eta_q_eps, link[m]); 
      }	     

      # Evaluate association structures
      if (has_assoc_ev[m] == 1) {
        mark = mark + 1;
	      e_eta_q = e_eta_q + a_beta[mark] * y_eta_q;
      }	
      if (has_assoc_cv[m] == 1) {
        mark = mark + 1;
        e_eta_q = e_eta_q + a_beta[mark] * ysep_q; 
      }	
      if (has_assoc_es[m] == 1) {
        vector[nrow_y_Xq] dydt_eta_q;
        dydt_eta_q = (y_eta_q_eps - y_eta_q) / eps;
        mark = mark + 1;
        e_eta_q = e_eta_q + a_beta[mark] * dydt_eta_q;          
      }
      if (has_assoc_cs[m] == 1) {
        vector[nrow_y_Xq] dydt_q;
        dydt_q = (ysep_q_eps - ysep_q) / eps;
        mark = mark + 1;
        e_eta_q = e_eta_q + a_beta[mark] * dydt_q;          
      }	
    }
  	if (sum_size_which_b > 0) {
  	  int mark_end;  // used to define segment of a_beta
  	  matrix[nrow_e_Xq,sum_size_which_b] x_assoc_sh;	  
  	  mark_end = mark + sum_size_which_b;
      mark = mark + 1;	  
  	  x_assoc_sh = make_x_assoc_sh(b, l, p, pmat, Npat, quadnodes, which_b,  
  	                               sum_size_which_b, size_which_b, t_i, M);
  	  e_eta_q = e_eta_q + x_assoc_sh * a_beta[mark:mark_end];
    }	
  }

  // Calculate log hazard at event times and unstandardised quadrature points 
  // NB assumes Weibull baseline hazard
  if (basehaz_weibull == 1) 
	  log_basehaz = log(weibull_shape[1]) + (weibull_shape[1] - 1) * e_log_times;
  else if (basehaz_splines == 1)
	  log_basehaz = e_ns_times * splines_coefs;	
  ll_haz_q = e_d .* (log_basehaz + e_eta_q);
					  
  // Partition event times and quad points
  ll_haz_eventtime = segment(ll_haz_q, 1, Npat);
  ll_haz_quadtime  = segment(ll_haz_q, (Npat + 1), Npat_times_quadnodes);
                
  // Log survival contribution to the likelihood (by summing over the 
  //   quadrature points to get the approximate integral) 
  ll_surv_eventtime = quadweight_times_half_eventtime .* exp(ll_haz_quadtime);        

  // Log likelihood for event model
  sum_ll_haz_eventtime = sum(ll_haz_eventtime);
  sum_ll_surv_eventtime = sum(ll_surv_eventtime);
  ll_event = sum_ll_haz_eventtime - sum_ll_surv_eventtime;				  

}
model {
  int disp_mark;
  int nois_mark;
  int int_mark;
  int int_markun;
  int int_marklo;
  int int_markup;
  vector[sum_y_N] y_eta;                                     

  disp_mark = 1;
  nois_mark = 1;
  int_mark = 1;
  int_markun = 1;
  int_marklo = 1;
  int_markup = 1;
  
  // Longitudinal submodel(s): regression equations
  if (sum_y_K > 0) y_eta = y_X * y_beta;
  else y_eta = rep_vector(0.0, sum_y_N);
  //if (y_has_offset == 1) y_eta = y_eta + y_offset;
  y_eta = y_eta + csr_matrix_times_vector(sum_y_N, len_b, w, v, u, b_by_model);
  for (m in 1:M) {
    if (y_has_intercept[m] == 1) {
      if (y_has_intercept_unbound[m] == 1)
        y_eta = y_eta + 
                    y_gamma_unbound[sum(y_has_intercept_unbound[1:m])];
      else if (y_has_intercept_lobound[m] == 1)
        y_eta = y_eta - min(y_eta) + 
                    y_gamma_lobound[sum(y_has_intercept_lobound[1:m])];
      else if (y_has_intercept_upbound[m] == 1)
        y_eta = y_eta - max(y_eta) + 
                    y_gamma_lobound[sum(y_has_intercept_upbound[1:m])];					
    }	
	
    if (y_centre == 1) {
      // correction to eta if model has no intercept (if X is centered)
      y_eta = y_eta + dot_product(y_xbar, y_beta); 
    }
    
#    if (family[m] == 8) {  # poisson-gamma mixture
#      if      (link[m] == 1) y_eta = y_eta + log(y_dispersion[disp_mark]) + log(y_noise[nois_mark]);
#      else if (link[m] == 2) y_eta = y_eta * y_dispersion[disp_mark] .* y_noise[nois_mark];
#      else                   y_eta = y_eta + sqrt(y_dispersion[disp_mark]) + sqrt_vec(y_noise[nois_mark]);
#    }    
    
    // Log-likelihood for longitudinal submodel(s)
    if (y_has_weights[m] == 0 && prior_PD == 0) { # unweighted log-likelihoods
      if (family[m] == 1) {
        if (link[m] == 1)      target += normal_lpdf(y_real | y_eta, y_dispersion[disp_mark]);
        else if (link[m] == 2) target += lognormal_lpdf(y_real | y_eta, y_dispersion[disp_mark]);
        else target += normal_lpdf(y_real | divide_real_by_vector(1, y_eta), y_dispersion[disp_mark]);
      }
      else if (family[m] == 2) {
        target += GammaReg(y_real, y_eta, y_dispersion[disp_mark], link[m], sum_log_y[m]);
      }
      else if (family[m] == 3) {
        target += inv_gaussian(y_real, 
                               linkinv_inv_gaussian(y_eta, link[m]), 
                               y_dispersion[disp_mark], sum_log_y[m], sqrt_y);
      }
	    else if (family[m] == 4) {
		    vector[y_N01[m,1]] y_eta0_tmp;
		    vector[y_N01[m,2]] y_eta1_tmp;
	      real dummy;  // irrelevant but useful for testing
		    y_eta0_tmp = segment(y_eta, 1, y_N01[m,1]);
		    y_eta1_tmp = segment(y_eta, (y_N01[m,1] + 1), y_N01[m,2]);
	      dummy = ll_bern_lp(y_eta0_tmp, y_eta1_tmp, link[m], y_N01[m,]);	  
	    }
	    else if (family[m] == 5) {
	      real dummy;  // irrelevant but useful for testing
        dummy = ll_binom_lp(y_int, trials, y_eta, link[m]);	  
	    }
	    else if (family[m] == 6 || family[m] == 8) {
        if (link[m] == 1) target += poisson_log_lpmf(y_int | y_eta);
        else target += poisson_lpmf(y_int | linkinv_count(y_eta, link[m]));
	    }
	    else if (family[m] == 7) {
  	    if (link[m] == 1) target += neg_binomial_2_log_lpmf(y_int | y_eta, y_dispersion[disp_mark]);
        else target += neg_binomial_2_lpmf(y_int | 
                                           linkinv_count(y_eta, link[m]), y_dispersion[disp_mark]);
	    }	    
    }    
    else if (prior_PD == 0) { # weighted log-likelihoods
  	  vector[y_N[m]] summands;
  	  if (family[m] == 1) {
  	    summands = pw_gauss(y_real, y_eta, y_dispersion[disp_mark], link[m]);
  	    target += dot_product(y_weights, summands);	  	    
  	  }
  	  else if (family[m] == 2) {
  	    summands = pw_gamma(y_real, y_eta, y_dispersion[disp_mark], link[m]);
  	    target += dot_product(y_weights, summands);	  	    
  	  }
  	  else if (family[m] == 3) {
  	    summands = pw_inv_gaussian(y_real, y_eta, y_dispersion[disp_mark], 
  		                           link[m], log_y, sqrt_y);
  	    target += dot_product(y_weights, summands);
  	  }
	    else if (family[m] == 4) {
    		vector[y_N01[m,1]] y_weights0_tmp;
    		vector[y_N01[m,2]] y_weights1_tmp;
    		vector[y_N01[m,1]] y_eta0_tmp;
    		vector[y_N01[m,2]] y_eta1_tmp;
    		y_eta0_tmp = segment(y_eta, 1, y_N01[m,1]);
    		y_eta1_tmp = segment(y_eta, (y_N01[m,1] + 1), y_N01[m,2]);
    		y_weights0_tmp = segment(y_weights, 1, y_N01[m,1]);
    		y_weights1_tmp = segment(y_weights, (y_N01[m,1] + 1), y_N01[m,2]);		
        target += dot_product(y_weights0_tmp, pw_bern(0, y_eta0_tmp, link[m]));
        target += dot_product(y_weights1_tmp, pw_bern(1, y_eta1_tmp, link[m]));
  	  }
  	  else if (family[m] == 5) {
        target += dot_product(y_weights, pw_binom(y_int, trials, y_eta, link[m]));
  	  }
  	  else if (family[m] == 6 || family[m] == 8) {
        target += dot_product(y_weights, pw_pois(y_int, y_eta, link[m]));    		
  	  }
  	  else if (family[m] == 7) {
        target += dot_product(y_weights, pw_nb(y_int, y_eta, y_dispersion[disp_mark], link[m]));
  	  }  	  
    }
    
    if (y_has_dispersion[m] == 1) disp_mark = disp_mark + 1;
#    if (y_has_noise[m] == 1)      nois_mark = nois_mark + 1;
  }  
                           
  // Log-likelihood for event submodel   
  if (prior_PD == 0) target += ll_event; 

  // Log-priors for coefficients in longitudinal submodel(s)
  if (priorLong_dist == 1) y_z_beta ~ normal(0, 1);
  else if (priorLong_dist == 2) {
    if (y_t_all_124) y_z_beta ~ normal(0,1);
    else if (y_t_any_124) for (k in 1:sum_y_K) {
      if (priorLong_df[k] == 1 || priorLong_df[k] == 2 || priorLong_df[k] == 4)
        y_z_beta[k] ~ normal(0,1);
      else y_z_beta[k] ~ student_t(priorLong_df[k], 0, 1);
    }
    else y_z_beta ~ student_t(priorLong_df, 0, 1);
  }
  else if (priorLong_dist == 3) { // hs
    y_z_beta ~ normal(0,1);
    y_local[1] ~ normal(0,1);
    y_local[2] ~ inv_gamma(0.5 * priorLong_df, 0.5 * priorLong_df);
    y_global[1] ~ normal(0,1);
    y_global[2] ~ inv_gamma(0.5, 0.5);
  }
  else if (priorLong_dist == 4) { // hs+
    y_z_beta ~ normal(0,1);
    y_local[1] ~ normal(0,1);
    y_local[2] ~ inv_gamma(0.5 * priorLong_df, 0.5 * priorLong_df);
    y_local[3] ~ normal(0,1);
    // unorthodox useage of priorLong_scale as another df hyperparameter
    y_local[4] ~ inv_gamma(0.5 * priorLong_scale, 0.5 * priorLong_scale);
    y_global[1] ~ normal(0,1);
    y_global[2] ~ inv_gamma(0.5, 0.5);
  }
  // else priorLong_dist is 0 and nothing is added 
  
  // Log-prior for intercept in longitudinal submodel(s)
  for (m in 1:M) {
    if (y_has_intercept[m] == 1) {
      if (y_has_intercept_unbound[m] == 1) { # unbound intercept
        if (priorLong_dist_for_intercept == 1)  // normal
          y_gamma_unbound[int_markun] ~ normal(priorLong_mean_for_intercept[int_mark], 
                                      priorLong_scale_for_intercept[int_mark]);
        else if (priorLong_dist_for_intercept == 2)  // student_t
          y_gamma_unbound[int_markun] ~ student_t(priorLong_df_for_intercept[int_mark], 
                                      priorLong_mean_for_intercept[int_mark], 
                                      priorLong_scale_for_intercept[int_mark]);
        // else priorLong_dist is 0 and nothing is added 
        int_markun = int_markun + 1;
      }
      else if (y_has_intercept_lobound[m] == 1)  { # lower bounded intercept
        if (priorLong_dist_for_intercept == 1)  // normal
          y_gamma_lobound[int_marklo] ~ normal(priorLong_mean_for_intercept[int_mark], 
                                    priorLong_scale_for_intercept[int_mark]);
        else if (priorLong_dist_for_intercept == 2)  // student_t
          y_gamma_lobound[int_marklo] ~ student_t(priorLong_df_for_intercept[int_mark], 
                                    priorLong_mean_for_intercept[int_mark], 
                                    priorLong_scale_for_intercept[int_mark]);
        // else priorLong_dist is 0 and nothing is added
        int_marklo = int_marklo + 1;
      }
      else if (y_has_intercept_lobound[m] == 1)  { # lower bounded intercept
        if (priorLong_dist_for_intercept == 1)  // normal
          y_gamma_upbound[int_markup] ~ normal(priorLong_mean_for_intercept[int_mark], 
                                    priorLong_scale_for_intercept[int_mark]);
        else if (priorLong_dist_for_intercept == 2)  // student_t
          y_gamma_upbound[int_markup] ~ student_t(priorLong_df_for_intercept[int_mark], 
                                    priorLong_mean_for_intercept[int_mark], 
                                    priorLong_scale_for_intercept[int_mark]);
        // else priorLong_dist is 0 and nothing is added
        int_markup = int_markup + 1;
      }
    int_mark = int_mark + 1;  
    }
  }
      
  // Log-priors for coefficients in event submodel
  if (priorEvent_dist == 1) e_z_beta ~ normal(0, 1);
  else if (priorEvent_dist == 2) {
    if (e_t_all_124) e_z_beta ~ normal(0,1);
    else if (e_t_any_124) for (k in 1:e_K) {
      if (priorEvent_df[k] == 1 || priorEvent_df[k] == 2 || priorEvent_df[k] == 4)
        e_z_beta[k] ~ normal(0,1);
      else e_z_beta[k] ~ student_t(priorEvent_df[k], 0, 1);
    }
    else e_z_beta ~ student_t(priorEvent_df, 0, 1);
  }
  else if (priorEvent_dist == 3) { // hs
    e_z_beta ~ normal(0,1);
    e_local[1] ~ normal(0,1);
    e_local[2] ~ inv_gamma(0.5 * priorEvent_df, 0.5 * priorEvent_df);
    e_global[1] ~ normal(0,1);
    e_global[2] ~ inv_gamma(0.5, 0.5);
  }
  else if (priorEvent_dist == 4) { // hs+
    e_z_beta ~ normal(0,1);
    e_local[1] ~ normal(0,1);
    e_local[2] ~ inv_gamma(0.5 * priorEvent_df, 0.5 * priorEvent_df);
    e_local[3] ~ normal(0,1);
    // unorthodox useage of priorEvent_scale as another df hyperparameter
    e_local[4] ~ inv_gamma(0.5 * priorEvent_scale, 0.5 * priorEvent_scale);
    e_global[1] ~ normal(0,1);
    e_global[2] ~ inv_gamma(0.5, 0.5);
  }
  // else priorEvent_dist is 0 and nothing is added 
  
  // Log-prior for intercept in event submodel
  if (e_has_intercept == 1) {
    if (priorEvent_dist_for_intercept == 1)  // normal
      e_gamma ~ normal(priorEvent_mean_for_intercept, priorEvent_scale_for_intercept);
    else if (priorEvent_dist_for_intercept == 2)  // student_t
      e_gamma ~ student_t(priorEvent_df_for_intercept, priorEvent_mean_for_intercept, 
                        priorEvent_scale_for_intercept);
    // else priorEvent_dist is 0 and nothing is added 
  }    

  // Log-priors for association parameters
  if (priorAssoc_dist == 1) a_z_beta ~ normal(0, 1);
  else if (priorAssoc_dist == 2) {
    if (a_t_all_124) a_z_beta ~ normal(0,1);
    else if (a_t_any_124) for (k in 1:a_K) {
      if (priorAssoc_df[k] == 1 || priorAssoc_df[k] == 2 || priorAssoc_df[k] == 4)
        a_z_beta[k] ~ normal(0,1);
      else a_z_beta[k] ~ student_t(priorAssoc_df[k], 0, 1);
    }
    else a_z_beta ~ student_t(priorAssoc_df, 0, 1);
  }
  else if (priorAssoc_dist == 3) { // hs
    a_z_beta ~ normal(0,1);
    a_local[1] ~ normal(0,1);
    a_local[2] ~ inv_gamma(0.5 * priorAssoc_df, 0.5 * priorAssoc_df);
    a_global[1] ~ normal(0,1);
    a_global[2] ~ inv_gamma(0.5, 0.5);
  }
  else if (priorAssoc_dist == 4) { // hs+
    a_z_beta ~ normal(0,1);
    a_local[1] ~ normal(0,1);
    a_local[2] ~ inv_gamma(0.5 * priorAssoc_df, 0.5 * priorAssoc_df);
    a_local[3] ~ normal(0,1);
    // unorthodox useage of priorAssoc_scale as another df hyperparameter
    a_local[4] ~ inv_gamma(0.5 * priorAssoc_scale, 0.5 * priorAssoc_scale);
    a_global[1] ~ normal(0,1);
    a_global[2] ~ inv_gamma(0.5, 0.5);
  }
  // else priorAssoc_dist is 0 and nothing is added   
    
  // Log-prior for scale(s)
  disp_mark = 1;
  for (m in 1:M) {
    if (y_has_dispersion[m] == 1) {
      if (priorLong_scale_for_dispersion[disp_mark] > 0) {
        target += cauchy_lpdf(y_dispersion_unscaled[disp_mark] | 0, 1);    
      }
    disp_mark = disp_mark + 1;  
    }
  }
   
  // Log-prior for Weibull shape
  if (basehaz_weibull == 1) 
    target += cauchy_lpdf(weibull_shape_unscaled | 0, 1);   
    
  // Log-prior for baseline hazard spline coefficients 
  if (basehaz_splines == 1) 
    target += normal_lpdf(splines_coefs_unscaled | 0, 1);
    
  // Prior for random effects model
  if (t > 0) decov_lp(z_b, z_T, rho, zeta, tau, 
                      regularization, delta, shape, t, p);
   
}
/*
generated quantities {
  real alpha[sum_y_has_intercept];
  real mean_PPD;
  mean_PPD = 0;
  if (has_intercept == 1)
    alpha[1] = gamma[1] - dot_product(xbar, beta);
  {
  vector[N] eta;  // linear predictor
  if (K > 0) eta = X * beta;
  else eta = rep_vector(0.0, N);
  if (has_offset == 1) eta = eta + offset;
    if (t > 0) eta = eta + csr_matrix_times_vector(N, len_b, w, v, u, b);
    if (has_intercept == 1) {
      if (family == 1 || link == 2) eta = eta + gamma[1];
      else {
        real min_eta;
        min_eta = min(eta);
        alpha[1] = alpha[1] - min_eta;
        eta = eta - min_eta + gamma[1];
      }
    }
    else {
      // correction to eta if model has no intercept (because X is centered)
      eta = eta + dot_product(xbar, beta); 
    }
    
    eta = segment(eta, 1, N);
    
    if (family == 1) {
      if (link > 1) eta = linkinv_gauss(eta, link);
      for (n in 1:N) mean_PPD = mean_PPD + normal_rng(eta[n], dispersion);
    }
    else if (family == 2) {
      if (link > 1) eta = linkinv_gamma(eta, link);
      for (n in 1:N) mean_PPD = mean_PPD + gamma_rng(dispersion, dispersion / eta[n]);
    }
    else if (family == 3) {
      if (link > 1) eta = linkinv_inv_gaussian(eta, link);
      for (n in 1:N) mean_PPD = mean_PPD + inv_gaussian_rng(eta[n], dispersion);
    }
    mean_PPD = mean_PPD / N;
  }
}
*/
