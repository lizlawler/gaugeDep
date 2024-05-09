#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec gauss_gauge(arma::vec const& W, double const& dep) {
  arma::vec w2 = 1 - W;
  return (W + w2 - 2 * dep * sqrt(W % w2)) / (1 - pow(dep, 2));
}

// [[Rcpp::export]]
double trunc_log_lhood(arma::vec const& W, arma::vec const& R, arma::vec const& r0_w,
                       double const& alpha, double const& theta) {
  int n = R.size();
  double log_lik = 0.0;
  
  arma::vec gw = gauss_gauge(W, theta);
  
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
arma::mat block_covariance(arma::mat param_block, arma::rowvec global_mean) {
  int n = param_block.n_rows;
  param_block.each_row() -= global_mean;
  // arma::mat cov = (1 / ((double)n - 1)) * X * X.t();
  arma::mat cov = param_block.t() * param_block;
  return cov/((double)n - 1);
}

// // [[Rcpp::export]]
// arma::mat trunc_am_lap(int const& n_iter, arma::vec const& W, arma::vec const& R,
//                         arma::vec const& r0_w,
//                         int const& dim=2,
//                         int const& steps=20,
//                         double const& r_opt=0.234) {
//   
//   // Initialize parameters from the priors
//   double alpha = R::rgamma(4.0, 0.5);
//   double theta = R::runif(0.0, 1.0);
//   
//   // Initialize memory storage for samples if required
//   arma::mat samples(n_iter, 4);
//   
//   // Initialize scale and covariance for proposal
//   double sigma_m = pow(2.38, 2) / (double)dim;
//   arma::mat cov = arma::eye(dim, dim);
//   
//   // Initialize block update parameters
//   int jumps = 0;
//   double c1 = 0.8;
//   double c0 = 1.0;
//   
//   arma::vec params_curr = {alpha, theta}; 
//   
//   // Sample for n iterations
//   for(int iter = 0; iter < n_iter; iter++) {
//     
//     arma::vec params_star = arma::mvnrnd(params_curr, sigma_m * cov);
//     double alpha_star = params_star(0);
//     double theta_star = params_star(1);
//     
//     double log_ratio = 0.0;
//     log_ratio += (log_prior_alpha(alpha_star) + log_prior_theta(theta_star) + trunc_log_lhood(W, R, r0_w, alpha_star, theta_star));
//     log_ratio -= (log_prior_alpha(alpha) + log_prior_theta(theta) + trunc_log_lhood(W, R, r0_w, alpha, theta));
//     
//     if(log(R::runif(0.0,1.0)) < log_ratio) {
//       alpha = alpha_star;
//       theta = theta_star;
//       jumps++;
//     }
//     
//     samples(iter, 0) = iter + 1;
//     samples(iter, 1) = alpha;
//     samples(iter, 2) = theta;
//     samples(iter, 3) = sigma_m;
//     
// 
//     if((iter + 1) % steps == 0 && iter != (n_iter - 1)) {
//       // Calculate block acceptance rate
//       double r_hat = jumps / (double)steps;
//       std::cout << "r_hat: " << r_hat << std::endl;
//       
//       // Subset samples to be the past 20 iterations
//       int end_row = iter;
//       int start_row = std::max(0, end_row - steps + 1);
//       arma::mat params_block = samples(arma::span(start_row, end_row), arma::span(1, 2));
// 
//       // Block covariance estimate 
//       arma::rowvec col_means = arma::mean(samples.head_rows(iter + 1).cols(1, 2), 0);
//       arma::mat cov_hat = block_covariance(params_block, col_means);
//       
//       // Adjust scaling and covariance
//       double block_num = (iter + 1) / (double)steps;
//       double gamma1 = pow(block_num + 5, -c1);
//       double gamma2 = c0 * gamma1;
//       sigma_m *= std::exp(gamma2 * (r_hat - r_opt));
//       cov *= (1 - gamma1);
//       cov += gamma1 * cov_hat;
//       
//       // Reset block acceptances 
//       jumps = 0;
//     }
//     
//     if((iter + 1) % 1000 == 0) {
//       // Check for user interrupt
//       Rcpp::checkUserInterrupt();
//     }
//   }
//   
//   return samples;
// }

// [[Rcpp::export]]
arma::mat trunc_am_lap_scaling(int const& n_iter, arma::vec const& W, arma::vec const& R,
                               arma::vec const& r0_w, arma::mat const& pilot_cov,
                               int const& dim=2,
                               int const& steps=20,
                               double const& r_opt=0.234) {
  
  // Initialize parameters from the priors
  double alpha = R::rgamma(4.0, 0.5);
  double theta = R::runif(0.0, 1.0);
  
  // Initialize memory storage for samples if required
  arma::mat samples(n_iter, 4);
  
  // Initialize scale and covariance for proposal
  double sigma_m = pow(2.38, 2) / (double)dim;
  
  // Initialize block update parameters
  int jumps = 0;
  double c1 = 0.8;
  double c0 = 1.0;
  
  arma::vec params_curr = {alpha, theta}; 
  
  // Sample for n iterations
  for(int iter = 0; iter < n_iter; iter++) {
    
    arma::vec params_star = arma::mvnrnd(params_curr, sigma_m * pilot_cov);
    double alpha_star = params_star(0);
    double theta_star = params_star(1);
    
    double log_ratio = 0.0;
    log_ratio += (log_prior_alpha(alpha_star) + log_prior_theta(theta_star) + trunc_log_lhood(W, R, r0_w, alpha_star, theta_star));
    log_ratio -= (log_prior_alpha(alpha) + log_prior_theta(theta) + trunc_log_lhood(W, R, r0_w, alpha, theta));
    
    if(log(R::runif(0.0,1.0)) < log_ratio) {
      alpha = alpha_star;
      theta = theta_star;
      jumps++;
    }
    
    samples(iter, 0) = iter + 1;
    samples(iter, 1) = alpha;
    samples(iter, 2) = theta;
    samples(iter, 3) = sigma_m;
    
    
    if((iter + 1) % steps == 0 && iter != (n_iter - 1)) {
      // Calculate block acceptance rate
      double r_hat = jumps / (double)steps;
      std::cout << "r_hat: " << r_hat << std::endl;
      
      // Adjust scaling and covariance
      double block_num = (iter + 1) / (double)steps;
      double gamma1 = pow(block_num + 5, -c1);
      double gamma2 = c0 * gamma1;
      sigma_m *= std::exp(gamma2 * (r_hat - r_opt));
      
      // Reset block acceptances 
      jumps = 0;
    }
    
    if((iter + 1) % 1000 == 0) {
      // Check for user interrupt
      Rcpp::checkUserInterrupt();
    }
  }
  
  return samples;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# library(RcppSimdJson)
# idx <- fload("data/gauss/low_11.json")$idx
# W <- fload("data/gauss/low_11.json")$W[idx]
# R <- fload("data/gauss/low_11.json")$R[idx]
# r0_w <- fload("data/gauss/low_11.json")$r0_w[idx]

samples <- trunc_mcmc_mh(10000, W, R, r0_w)
pilot_cov <- cov(samples[,2:3])
samples_am <- trunc_am_lap_scaling(10000, W, R, r0_w, pilot_cov = pilot_cov, steps = 20)
# sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = TRUE)
# X <- mvtnorm::rmvnorm(10, c(0.5,1.5), sigma = sigma)
# block_covariance(X, colMeans(X[1:5,]))
# custom_covariance(t(X))

*/
