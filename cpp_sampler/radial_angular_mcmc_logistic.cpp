#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::vec gauss_gauge(arma::vec const& W, double const& dep) {
  vec w2 = 1 - W;
  return (W + w2 - 2 * dep * sqrt(W % w2)) / (1 - pow(dep, 2));
}

// [[Rcpp::export]]
arma::vec logistic_gauge(arma::vec const& W, double const& dep) {
  mat w1_and_2(W.size(), 2);
  w1_and_2.col(0) = W;
  w1_and_2.col(1) = 1-W;
  double r_inv = 1/dep;
  return(r_inv * max(w1_and_2, 1) + (1 - r_inv) * min(w1_and_2, 1));
}

// [[Rcpp::export]]
double est_vol(arma::vec const& sum_term, arma::vec const& sqrt_term, double const& pars) {
  vec gx = (sum_term - 2 * pars * sqrt_term) / (double)(1 - pow(pars, 2));
  return mean(conv_to<vec>::from(gx <= 1));
}

// [[Rcpp::export]]
double angular_loglik(arma::vec const& W, double const& pars,
                      arma::vec const& sum_term, arma::vec const& sqrt_term, int const& dim) {
  int N = W.size();
  double L_volume = est_vol(sum_term, sqrt_term, pars);
  return(-(double)dim * sum(log(gauss_gauge(W, pars))) - (double)N * (log((double)dim) + log(L_volume)));
}

// [[Rcpp::export]]
double radial_cens_loglik(arma::vec const& radii, arma::vec const& pars,
                          arma::vec const& threshold, arma::vec const& W) {
  int n = radii.size();
  double alpha = pars(0);
  double dep = pars(1);
  vec beta = logistic_gauge(W, dep);
  double loglik = 0.0;
  for(int i = 0; i < n; i ++) {
    if(radii(i) < threshold(i)) {
      loglik += R::pgamma(threshold(i), alpha, 1/beta(i), true, true);
    } else {
      loglik += R::dgamma(radii(i), alpha, 1/beta(i), true);
    }
  }
  return(loglik);
}

// [[Rcpp::export]]
double joint_loglik(Rcpp::List const& data_list, arma::vec const& params, arma::vec const& threshold,
                    arma::vec const& sum_term, arma::vec const& sqrt_term, int const& dim) {
  vec R = data_list["R"];
  vec W = data_list["W"];
  double angular_ll = angular_loglik(W, params(1), sum_term, sqrt_term, dim);
  double radial_ll = radial_cens_loglik(R, {params(0), params(2)}, threshold, W);
  return angular_ll + radial_ll;
}

// [[Rcpp::export]]
double prior_fn(arma::vec const& params) {
  // double alpha_prior = R::dgamma(params(0), hyper_params(0), 1/hyper_params(1), true);
  double alpha_prior = R::dgamma(params(0), 4, 0.5, true);
  double ang_dep_prior = R::dunif(params(1), 0, 1, true);
  double rad_dep_prior = R::dunif(params(2), 0, 1, true);
  return alpha_prior + ang_dep_prior + rad_dep_prior;
}

// [[Rcpp::export]]
Rcpp::List adaptive_mh_rad_ang_vol(Rcpp::List const& z,
                                   arma::vec const& threshold,
                                   arma::vec const& sum_term,
                                   arma::vec const& sqrt_term,
                                   int const& dim,
                                   arma::rowvec const& starting_theta,
                                   int const& n_updates,
                                   int const& update_freq,
                                   Rcpp::Nullable<Rcpp::NumericMatrix> prop_Sigma_ = R_NilValue,
                                   bool const& adapt_cov = false,
                                   double const& r_opt = 0.234,
                                   double const& c0 = 10,
                                   double const& c1 = 0.8,
                                   int const& K = 10) {

  double eps = 0.001;
  int k = 3;

  int p = starting_theta.size();

  mat prop_Sigma;
  // Check if proposal matrix is provided and valid
  if (prop_Sigma_.isNotNull()) {
    prop_Sigma = Rcpp::as<arma::mat>(prop_Sigma_);
    if (prop_Sigma.size() != p * p || !prop_Sigma.is_symmetric()) {
      Rcpp::Rcout << "Proposal matrix is invalid, using identity matrix." << std::endl;
      prop_Sigma = arma::eye(p, p);  // Use identity matrix
    }
  } else {
    Rcpp::Rcout << "Proposal matrix not provided, using identity matrix." << std::endl;
    prop_Sigma = eye(p, p);  // Use identity matrix
  }

  // Initialize sigma_m to rule of thumb
  double sigma_m = pow(2.4 / (double)p, 2);
  double r_hat = 0;

  // Initialize cholesky decomp of proposal covariance matrix
  mat prop_C = chol(prop_Sigma);

  // Setup and initialize trace/sample objects
  mat trace = mat(n_updates, p, fill::zeros);
  vec sigma_m_trace = vec(n_updates, fill::zeros);
  vec r_trace = vec(n_updates, fill::zeros);
  vec jump_trace = vec(n_updates, fill::zeros);

  trace.row(0) = starting_theta;
  sigma_m_trace(0) = sigma_m;

  // Initialize MH
  vec theta = starting_theta.t();
  double likelihood = joint_loglik(z, theta, threshold, sum_term, sqrt_term, dim);
  double prior = prior_fn(theta);

  ///////////////////////////////////////////////////////////////////////////////
  // Begin main loop
  for(int iter = 1; iter < n_updates; iter++) {
    vec theta_star = theta + sigma_m * (randn<rowvec>(p) * prop_C).t();
    double prior_star = prior_fn(theta_star);
    if(prior_star != -datum::inf) {
      double likelihood_star = joint_loglik(z, theta_star, threshold, sum_term, sqrt_term, dim);
      double mh_ratio = exp(prior_star + likelihood_star - prior - likelihood);
      if(mh_ratio > randu()) {
        theta = theta_star;
        prior = prior_star;
        likelihood = likelihood_star;
        jump_trace(iter) = 1;
      }
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Adaptive piece
    if(((iter + 1) % K) == 0) {
      int iter_star = iter + 1;
      double gamma1 = c0 / pow(((iter_star / K) + k), c1);
      double gamma2 = 1 / pow(((iter_star / K) + k), c1);

      r_hat = mean(jump_trace.subvec((iter_star - K), iter));
      sigma_m = exp(log(sigma_m) + gamma1 * (r_hat - r_opt));

      if(adapt_cov) {
        prop_Sigma = prop_Sigma +
          gamma2 * (cov(trace(span(iter_star - K, iter), span(0, (p - 1)))) - prop_Sigma);

        while(!prop_Sigma.is_symmetric()) {
          prop_Sigma = prop_Sigma + eps * eye(p, p);
          Rcpp::Rcout << "Oops, proposal covariance matrix is now:\n" << std::endl;
          Rcpp::Rcout << prop_Sigma << std::endl;
        }
        prop_C = chol(prop_Sigma);
      }
    }

    // End adaptation phase                                                      //
    ///////////////////////////////////////////////////////////////////////////////

    // Update trace/sample objects
    trace.row(iter) = theta.t();
    sigma_m_trace(iter) = sigma_m;
    r_trace(iter) = r_hat;

    // Report progress to console
    if(((iter + 1) % update_freq) == 0) {
      Rcpp::Rcout << "Finished " << (iter + 1) << " out of " << n_updates << " iterations." << std::endl;
      Rcpp::checkUserInterrupt();
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("trace") = trace,
    Rcpp::Named("sigma_m_trace") = sigma_m_trace,
    Rcpp::Named("r_trace") = r_trace,
    Rcpp::Named("acc_prob") = mean(jump_trace)
  );
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
data <- RcppSimdJson::fload("data/logistic/low_1.json")
W <- data$W
R <- data$R
r0w <- data$r0_w
data_list <- list(R = R, W = W)
grid <- expand.grid(seq(0, 1, length.out=100), seq(0, 1, length.out=100))
sum_grid <- grid[,1] + grid[,2]
sqrt_grid <- sqrt(grid[,1] * grid[,2])
samples <- adaptive_mh_rad_ang_vol(z = data_list,
                                 threshold = r0w,
                                 sum_term = sum_grid,
                                 sqrt_term = sqrt_grid,
                                 dim = 2,
                                 starting_theta = c(rgamma(1, 4, 2), runif(2)),
                                 n_updates = 25000, update_freq = 500,
                                 adapt_cov = TRUE,
                                 K = 10)
*/
