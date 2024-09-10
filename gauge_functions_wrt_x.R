## gauge_wrt_x functions
gauss_gauge_wrt_x <- function(x1, x2, dep_par = 0.5) {
  top <- x1 + x2 - 2 * dep_par * sqrt(x1 * x2)
  return(top/(1-dep_par^2))
}

logistic_gauge_wrt_x <- function(x1, x2, dep_par = 0.5) {
  r_inv <- 1/dep_par
  return(r_inv * pmax(x1, x2) + (1-r_inv)*pmin(x1, x2))
}

inv_log_gauge_wrt_x <- function(x1, x2, dep_par = 0.5) ((x1^(1/dep_par) + x2^(1/dep_par))^dep_par)

asym_log_gauge_wrt_x <- function(x1, x2, dep_par = 0.5) {
  r_inv <- 1/dep_par
  return(pmin((x1 + x2), (r_inv * pmax(x1, x2) + (1-r_inv) * pmin(x1, x2))))
}

dirichlet_gauge_wrt_x <- function(x1, x2, dep_par) {
  theta1 <- dep_par[1]
  theta2 <- dep_par[2]
  return((1 + theta1 + theta2) * pmax(x1, x2) - (theta1 * x1 + theta2 * x2))
}

rectangular_gauge_wrt_x <- function(x1, x2, dep_par) {
  return(pmax((x1 - x2) / dep_par, (x2 - x1) / dep_par, (x1 + x2) / (2 - dep_par)))
}