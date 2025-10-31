#' Computes the CIs and probabilities
#' @param x Numeric vector (endogenous variable).
#' @param y Numeric vector (outcome variable).
#' @param z Numeric vector (Instrumental variable).
#' @param alpha Numeric scalar for the level of H0 test.
#' @param rxu_range Numeric interval for r_xu values.
#' @return Simple CIs.
#' @export

ci_simple_union <- function(x, y, z,
                            rxu_range=c(0,0.8),
                            alpha = 0.05) {
  D_grid = seq(rxu_range[1], rxu_range[2], length.out = 80)
  n <- length(x)
  rho_xy <- cor(x,y,use = "complete.obs",method=c("pearson"))
  rho_xz <- cor(x,z,use = "complete.obs",method=c("pearson"))
  rho_yz <- cor(y,z,use = "complete.obs",method=c("pearson"))

  # Variances (simplified)
  var_xy <- (1 - rho_xy^2)^2 / n
  var_xz <- (1 - rho_xz^2)^2 / n
  var_yz <- (1 - rho_yz^2)^2 / n
  Sigma_rho <- diag(c(var_xy, var_xz, var_yz))

  g_vals <- sapply(D_grid, function(r) g_fun(r, rho_xy, rho_xz, rho_yz))
  idx_l <- which.min(g_vals); idx_u <- which.max(g_vals)
  plug_in <- c(g_vals[idx_l], g_vals[idx_u])

  # gradient at those r
  grad_l <- g_grad(D_grid[idx_l], rho_xy, rho_xz, rho_yz)
  grad_u <- g_grad(D_grid[idx_u], rho_xy, rho_xz, rho_yz)

  var_l <- t(grad_l) %*% Sigma_rho %*% grad_l
  var_u <- t(grad_u) %*% Sigma_rho %*% grad_u

  K <- length(D_grid)
  crit <- qnorm(1 - alpha / (2 * K * 2))  # conservative two-sided
  CI_lower <- plug_in[1] - crit * sqrt(var_l)
  CI_upper <- plug_in[2] + crit * sqrt(var_u)

  list(CI = c(CI_lower, CI_upper), plug_in = plug_in)
}
