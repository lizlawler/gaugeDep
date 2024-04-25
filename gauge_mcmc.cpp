#include <RcppArmadillo.h>
#include <boost/math/distributions/gamma.hpp>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
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
double cens_log_lhood(arma::vec const& W, arma::vec const& R, arma::vec const& r0_w,
                      double const& alpha, double const& theta) {
  int n = R.size();
  double log_lik = 0.0;
  
  arma::vec gw = gauss_gauge(W, theta);
  
  for(int i = 0; i < n; i++) {
    if(R(i) < r0_w(i)) {
      double log_cdf = R::pgamma(r0_w(i), alpha, 1/gw(i), true, true);
      log_lik += log_cdf;
    } else {
      double log_pdf = R::dgamma(R(i), alpha, 1/gw(i), true);
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
void trunc_mcmc_mh(int const& n_iter, arma::vec const& W, arma::vec const& R, arma::vec const& r0_w, 
                         arma::vec const& step_size, 
                         double const& update_int=10,
                         std::string const& output_file="output.csv") {

  auto start_time = std::chrono::high_resolution_clock::now();
  
  // Get the current date
  std::time_t now = std::time(nullptr);
  std::tm*local_time = std::localtime(&now);
  int year = local_time->tm_year + 1900;
  int month = local_time->tm_mon + 1;
  int day = local_time->tm_mday;
  int hour = local_time->tm_hour;
  int min = local_time->tm_min;
  int sec = local_time->tm_sec;

  // Initialize progress update variables
  int progress_update_interval = n_iter / update_int; // Update progress every (n_iter/update_int/100)% of iterations
  int next_progress_update = progress_update_interval;
  int progress = 0;
  
  // Set the width of the progress bar
  int bar_width = 70;
  
  // Initialize parameters from the priors
  double alpha = R::rgamma(4.0, 0.5);
  double theta = R::runif(0.0, 1.0);
  
  // Create storage for posterior log-likelihood
  double posterior_log_lik;
  
  // Open the output file for writing
  std::ofstream outfile(output_file);
  
  // Write metadata to the CSV file
  outfile << "metadata\n";
  outfile << "Date:" << year << "_" << month << "_" << day << "_" << hour << "_" << min << "_" << sec << "\n";
  outfile << "Initial values: alpha = " << alpha << ", theta = " << theta << "\n";
  outfile << "Sampling settings: num_iterations = " << n_iter << ", step_size = " << step_size << "...\n";
  outfile << "end_metadata\n";
  
  // Write column headers for the CSV file
  outfile << "iter,alpha,theta,posterior_log_lik\n";
  
  // Keep track of acceptances
  int accept = 0;
  
  // Sample for n iterations
  for(int iter = 0; iter < n_iter; ++iter) {
    
    // Check for user interrupt
    if (iter % 500 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    // Update progress
    if (iter + 1 == next_progress_update) {
      // Calculate percentage complete
      progress = (iter + 1) * 100 / n_iter;
      
      // Print progress bar
      std::cout << "[" << std::string(progress * bar_width / 100, '=') << std::string(bar_width - progress * bar_width / 100, ' ') << "] " << progress << "% (" << iter + 1 << "/" << n_iter << ")\r";
      std::cout.flush();
      
      // Update next progress update
      next_progress_update += progress_update_interval;
    }
    
    double alpha_star = alpha + R::rnorm(0.0, step_size(0));
    double theta_star = theta + R::rnorm(0.0, step_size(1));
    
    double log_ratio = 0.0;
    log_ratio += log_prior_alpha(alpha_star);
    log_ratio += log_prior_theta(theta_star);
    log_ratio += trunc_log_lhood(W, R, r0_w, alpha_star, theta_star);
    log_ratio -= log_prior_alpha(alpha);
    log_ratio -= log_prior_theta(theta);
    log_ratio -= trunc_log_lhood(W, R, r0_w, alpha, theta);
    
    if(log(R::runif(0.0,1.0)) < log_ratio) {
      alpha = alpha_star;
      theta = theta_star;
      accept++;
    }
    
    posterior_log_lik = trunc_log_lhood(W, R, r0_w, alpha, theta);
    
    // Write output to csv at each iteration
    // outfile << (iter + 1) << "," << alpha_samples(iter) << "," << theta_samples(iter) << "," << posterior_log_lik(iter) << "\n";
    outfile << (iter + 1) << "," << alpha << "," << theta << "," << posterior_log_lik << "\n";
  }

  // Print newline after the progress bar is complete
  std::cout << std::endl;
  
  // Compute elapsed time
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end_time - start_time;
  double total_time = elapsed_seconds.count();
  
  double accept_rate = accept / (double)n_iter;
  
  // Write elapsed time and acceptance rate to the output file
  outfile << "metadata\n";
  outfile << "Elapsed time: " << total_time << " seconds\n";
  outfile << "Acceptance rate: " << accept_rate << "\n";
  outfile << "end_metadata\n";
  
  // Close the output file
  outfile.close();
}

// [[Rcpp::export]]
Rcpp::List cens_mcmc_mh(int const& n_iter, arma::vec const& W, arma::vec const& R, arma::vec const& r0_w, 
                        arma::vec const& step_size, 
                        double const& update_int=10,
                        std::string const& output_file="output.csv") {
  auto start_time = std::chrono::high_resolution_clock::now();
  
  // Get the current date
  std::time_t now = std::time(nullptr);
  std::tm*local_time = std::localtime(&now);
  int year = local_time->tm_year + 1900;
  int month = local_time->tm_mon + 1;
  int day = local_time->tm_mday;
  int hour = local_time->tm_hour;
  int min = local_time->tm_min;
  int sec = local_time->tm_sec;
  
  // Initialize progress update variables
  int progress_update_interval = n_iter / update_int; // Update progress every (n_iter/update_int/100)% of iterations; defaults to 10%
  int next_progress_update = progress_update_interval;
  int progress = 0;
  
  // Set the width of the progress bar
  int bar_width = 70;
  
  // Initialize parameters from the priors
  double alpha = R::rgamma(4.0, 0.5);
  double theta = arma::randu(arma::distr_param(0.0, 1.0));
  
  // Create storage for samples
  arma::vec alpha_samples(n_iter);
  arma::vec theta_samples(n_iter);
  arma::vec posterior_log_lik(n_iter);
  
  // Keep track of acceptances
  int accept = 0;
  
  // Sample for n iterations
  for(int iter = 0; iter < n_iter; ++iter) {
    
    // Check for user interrupt every 500 iterations
    if (iter % 500 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    // Update progress
    if (iter + 1 == next_progress_update) {
      // Calculate percentage complete
      progress = (iter + 1) * 100 / n_iter;
      
      // Print progress bar
      std::cout << "[" << std::string(progress * bar_width / 100, '=') << std::string(bar_width - progress * bar_width / 100, ' ') << "] " << progress << "% (" << iter + 1 << "/" << n_iter << ")\r";
      std::cout.flush();
      
      // Update next progress update
      next_progress_update += progress_update_interval;
    }
    
    double alpha_star = alpha + R::rnorm(0.0, step_size(0));
    double theta_star = theta + R::rnorm(0.0, step_size(1));
    
    double log_ratio = 0.0;
    log_ratio += log_prior_alpha(alpha_star);
    log_ratio += log_prior_theta(theta_star);
    log_ratio += cens_log_lhood(W, R, r0_w, alpha_star, theta_star);
    log_ratio -= log_prior_alpha(alpha);
    log_ratio -= log_prior_theta(theta);
    log_ratio -= cens_log_lhood(W, R, r0_w, alpha, theta);
    
    if(log(R::runif(0.0,1.0)) < log_ratio) {
      alpha = alpha_star;
      theta = theta_star;
      accept++;
    }
    
    alpha_samples(iter) = alpha;
    theta_samples(iter) = theta;
    posterior_log_lik(iter) = cens_log_lhood(W, R, r0_w, alpha, theta);
  }
  
  // Print newline after the progress bar is complete
  std::cout << std::endl;
  
  // Compute elapsed time
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end_time - start_time;
  double total_time = elapsed_seconds.count();
  
  double accept_rate = accept / (double)n_iter;
  
  return List::create(_["alpha_samples"] = alpha_samples,
                      _["theta_samples"] = theta_samples,
                      _["posterior_log_lik"] = posterior_log_lik,
                      _["accept_rate"] = accept_rate,
                      _["total_time"] = total_time);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
idx <- RcppSimdJson::fload("data/gauss/low_10.json")$idx
W <- RcppSimdJson::fload("data/gauss/low_10.json")$W
R <- RcppSimdJson::fload("data/gauss/low_10.json")$R
r0_w <- RcppSimdJson::fload("data/gauss/low_10.json")$r0_w
# 
# # test <- trunc_log_lhood(W, R, r0_w, 2, 0.1)
trunc_mcmc_mh(5000, W[idx], R[idx], r0_w[idx], c(0.1, 0.2), update_int = 20, output_file="./test.csv")
# test_run <- cens_mcmc_mh(5000, W, R, r0_w, c(0.05, 0.05), update_int = 20)
*/
