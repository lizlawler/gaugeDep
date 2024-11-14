## Gauge functions
gauss_gauge <- function(w, dep_par = 0.5) {
  if(!is.null(nrow(w))) {
    d <- ncol(w)
    sigma <- matrix(dep_par, d, d)
    diag(sigma) <- 1
    sigma_inv <- solve(sigma)
    sqrt_w <- sqrt(w)
    return(rowSums((sqrt_w %*% sigma_inv) * sqrt_w))
  } else {
    top <- w + (1 - w) - 2 * dep_par * sqrt(w * (1 - w))
    return(top/(1-dep_par^2))
  }
}

w1 <- seq(0, 1, length.out = 200)
gw <- gauss_gauge(w1, 0.9)
plot(w1 / gw, (1-w1)/gw, pch = 20)

w1 <- runif(7)
w2 <- runif(7)
w3 <- runif(7)
w_mat <- cbind(w1, w2)
row_sums_w <- rowSums(w_mat)
w1 <- w1 / row_sums_w
w2 <- w2 / row_sums_w
w3 <- w3 / row_sums_w
w_mat <- cbind(w1, w2) |> as.matrix(nrow = 15, ncol = 2)
gauss_gauge(w_mat, 0.1)
min_vals <- apply(w_mat, 1, min)

logistic_gauge <- function(w, dep_par = 0.5) {
  gamma_inv <- 1 / dep_par
  if(!is.null(nrow(w))) {
    d <- ncol(w)
    min_vec <- apply(w, 1, min)
    return(gamma_inv * rowSums(w_mat) + (1 - d * gamma_inv) * min_vec)
  } else {
    return(gamma_inv * pmax(w, (1 - w)) + (1-gamma_inv)*pmin(w,(1 - w)))
  }
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