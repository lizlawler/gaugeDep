#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec gauss_gauge(arma::vec const& w1, double const& dep) {
  arma::vec w2 = 1 - w1;
  return (w1 + w2 - 2 * dep * sqrt(w1 % w2)) / (1 - pow(dep, 2));
}

// [[Rcpp::export]]
arma::vec logistic_gauge(arma::vec const& w1, double const& dep) {
  arma::vec w2 = 1 - w1;
  double dep_inv = 1/dep;
  return (dep_inv * arma::max(w1, w2) + (1-dep_inv) * arma::min(w1,w2));
}

// [[Rcpp::export]]
double trunc_log_lhood(arma::vec const& W, arma::vec const& R, arma::vec const& r0_w,
                       double const& alpha, double const& theta, bool gauss=true) {
  int n = R.size();
  double log_lik = 0.0;
  arma::vec gw(n);
  
  if(gauss) {
    gw = gauss_gauge(W, theta); 
  } else{
    gw = logistic_gauge(W, theta);
  }
  
  for(int i = 0; i < n; i++) {
    double gw_inv = 1/gw(i);
    double log_pdf = R::dgamma(R(i), alpha, gw_inv, true);
    double log_ccdf = R::pgamma(r0_w(i), alpha, gw_inv, false, true);
    log_lik += log_pdf;
    log_lik -= log_ccdf;
  }
  return log_lik;
}

// [[Rcpp::export]]
double cens_log_lhood(arma::vec const& W, arma::vec const& R, arma::vec const& r0_w,
                      double const& alpha, double const& theta, bool gauss=true) {
  int n = R.size();
  double log_lik = 0.0;
  arma::vec gw(n);
  
  if(gauss) {
    gw = gauss_gauge(W, theta); 
  } else{
    gw = logistic_gauge(W, theta);
  }

  for(int i = 0; i < n; i++) {
    double gw_inv = 1/gw(i);
    if(R(i) < r0_w(i)) {
      double log_cdf = R::pgamma(r0_w(i), alpha, gw_inv, true, true);
      log_lik += log_cdf;
    } else {
      double log_pdf = R::dgamma(R(i), alpha, gw_inv, true);
      log_lik += log_pdf;
    }
  }
  return log_lik;
}

// [[Rcpp::export]]
double log_prior_theta(double const& theta) {
  if(0 <= theta && theta <= 1) {
    return 0.0;
  } else{
    return R_NegInf;
  }
}

// [[Rcpp::export]]
double log_prior_alpha(double const& alpha) {
  return R::dgamma(alpha, 4, 0.5, true);
}

// [[Rcpp::export]]
arma::mat mcmc_mh(int const& n_iter, arma::vec const& W, arma::vec const& R,
                        arma::vec const& r0_w, arma::vec const& step_size, 
                        bool trunc=true, bool gauss=true) {

  // Initialize parameters from the priors
  double alpha = R::rgamma(4.0, 0.5);
  double theta = R::runif(0.0, 1.0);

  // Initialize memory storage for samples if required
  arma::mat samples(n_iter, 3);
  
  int accepts = 0;

  // Sample for n iterations
  for(int iter = 0; iter < n_iter; iter++) {
    
    double alpha_star = alpha + R::rnorm(0, step_size(0));
    double theta_star = theta + R::rnorm(0, step_size(1));

    double log_ratio = 0.0;
    log_ratio += log_prior_alpha(alpha_star) + log_prior_theta(theta_star);
    log_ratio -= log_prior_alpha(alpha) - log_prior_theta(theta);
    
    if(trunc) {
      if(gauss) {
        log_ratio += trunc_log_lhood(W, R, r0_w, alpha_star, theta_star, true);
        log_ratio -= trunc_log_lhood(W, R, r0_w, alpha, theta, true);
      } else {
        log_ratio += trunc_log_lhood(W, R, r0_w, alpha_star, theta_star, false);
        log_ratio -= trunc_log_lhood(W, R, r0_w, alpha, theta, false);
      } 
    } else {
      if(gauss) {
        log_ratio += cens_log_lhood(W, R, r0_w, alpha_star, theta_star, true);
        log_ratio -= cens_log_lhood(W, R, r0_w, alpha, theta, true);
      } else {
        log_ratio += cens_log_lhood(W, R, r0_w, alpha_star, theta_star, false);
        log_ratio -= cens_log_lhood(W, R, r0_w, alpha, theta, false);
      } 
    }

    if(log(R::runif(0.0,1.0)) < log_ratio) {
      alpha = alpha_star;
      theta = theta_star;
      accepts++;
    }

    samples(iter, 0) = iter + 1;
    samples(iter, 1) = alpha;
    samples(iter, 2) = theta;

    if((iter + 1) % 1000 == 0) {
      // Check for user interrupt
      Rcpp::checkUserInterrupt();
    }
  }
  std::cout << "Acceptance rate: " << accepts / (double)n_iter << std::endl;
  return samples;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
