## Gauge functions
gauss_gauge <- function(w, dep_par = 0.5) {
  top <- w + (1 - w) - 2 * dep_par * sqrt(w * (1 - w))
  return(top/(1-dep_par^2))
}

logistic_gauge <- function(w, dep_par = 0.5) {
  r_inv <- 1/dep_par
  return(r_inv * pmax(w, (1 - w)) + (1-r_inv)*pmin(w,(1 - w)))
}

inv_log_gauge <- function(w, dep_par = 0.5) ((w^(1/dep_par) + (1 - w)^(1/dep_par))^dep_par)

asym_log_gauge <- function(w, dep_par = 0.5) {
  r_inv <- 1/dep_par
  return(pmin((w + (1 - w)), (r_inv * pmax(w, (1 - w)) + (1-r_inv)*pmin(w,(1 - w)))))
}

dirichlet_gauge <- function(w, dep_par) {
  theta1 <- dep_par[1]
  theta2 <- dep_par[2]
  return((1 + theta1 + theta2) * pmax(w, (1 - w)) - (theta1 * w + theta2 * (1 - w)))
}

rectangular_gauge <- function(w, dep_par) {
  return(pmax((w - (1 - w)) / dep_par, ((1 - w) - w) / dep_par, (w + (1 - w)) / (2 - dep_par)))
}