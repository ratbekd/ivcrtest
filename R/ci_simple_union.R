#' Computes the CIs and probabilities
#' @param x Numeric vector (endogenous variable).
#' @param y Numeric vector (outcome variable).
#' @param z Numeric vector (Instrumental variable).
#' @param alpha Numeric scalar for the level of H0 test.
#' @param rxu_range Numeric interval for r_xu values.
#' @param grid_length Numeric scalar for the grid length.
#' @param cov_method Character value to distinguish between bootsrap and delta methods.
#' @param B_boot Numerical scalar for the bootsrap size.
#' @param seed Numerical scalar for the seed of the replications.
#' @param tol Numerical scalar for the tolerance level in bootstrap.
#' @param make_psd  Logical value.
#' @param Sigma_rho a covariance matrix.
#' @return Simple CIs.
#' @export

ci_simple_union <- function(x, y, z,
                            rxu_range   = c(0, 0.8),
                            grid_length = 80,
                            alpha       = 0.05,
                            # If you already computed a full 3x3 cov of (rho_xy,rho_xz,rho_yz), pass it here:
                            Sigma_rho   = NULL,
                            cov_method  = c("bootstrap", "diag"),
                            B_boot      = 800,
                            seed        = 123,
                            tol         = 1e-10,
                            make_psd    = TRUE) {

  cov_method <- match.arg(cov_method)

  # -----------------------------
  # 0) Clean data
  # -----------------------------
  dat <- data.frame(x = x, y = y, z = z)
  dat <- dat[complete.cases(dat), , drop = FALSE]
  n <- nrow(dat)
  if (n < 30) stop("ci_simple_union: need at least 30 complete observations.")

  # -----------------------------
  # 1) Sample correlations (rho_hat)
  # Order: (rho_xy, rho_xz, rho_yz)
  # -----------------------------
  rho_xy <- suppressWarnings(cor(dat$x, dat$y))
  rho_xz <- suppressWarnings(cor(dat$x, dat$z))
  rho_yz <- suppressWarnings(cor(dat$y, dat$z))
  deltahat <- c(rho_xy, rho_xz, rho_yz)
  names(deltahat) <- c("rho_xy", "rho_xz", "rho_yz")

  if (any(!is.finite(deltahat))) {
    stop("ci_simple_union: non-finite sample correlations (check for constant variables).")
  }

  # Clamp helper (avoid division by zero in sqrt(1 - rho^2))
  clamp <- function(a) pmin(pmax(a, -1 + tol), 1 - tol)
  rho_xy_c <- clamp(rho_xy)

  # -----------------------------
  # 2) Full covariance of correlation vector
  # -----------------------------
  if (is.null(Sigma_rho)) {
    if (cov_method == "bootstrap") {
      set.seed(seed)
      Rb <- matrix(NA_real_, nrow = B_boot, ncol = 3)
      for (b in seq_len(B_boot)) {
        idx <- sample.int(n, n, replace = TRUE)
        xb <- dat$x[idx]; yb <- dat$y[idx]; zb <- dat$z[idx]
        Rb[b, ] <- c(
          suppressWarnings(cor(xb, yb)),
          suppressWarnings(cor(xb, zb)),
          suppressWarnings(cor(yb, zb))
        )
      }
      Rb <- Rb[complete.cases(Rb), , drop = FALSE]
      if (nrow(Rb) < 50) stop("ci_simple_union: too many NA bootstrap draws; check data variation.")
      Sigma_rho <- cov(Rb)
      Sigma_rho <- (Sigma_rho + t(Sigma_rho)) / 2

      if (make_psd && requireNamespace("Matrix", quietly = TRUE)) {
        Sigma_rho <- as.matrix(Matrix::nearPD(Sigma_rho, corr = FALSE)$mat)
      }
    } else {
      # Fallback: diagonal-only (NOT consistent with "full covariance" claim)
      var_xy <- (1 - rho_xy^2)^2 / n
      var_xz <- (1 - rho_xz^2)^2 / n
      var_yz <- (1 - rho_yz^2)^2 / n
      Sigma_rho <- diag(c(var_xy, var_xz, var_yz))
    }
  } else {
    # sanity checks
    if (!is.matrix(Sigma_rho) || any(dim(Sigma_rho) != c(3, 3))) {
      stop("ci_simple_union: Sigma_rho must be a 3x3 matrix for (rho_xy, rho_xz, rho_yz).")
    }
    Sigma_rho <- (Sigma_rho + t(Sigma_rho)) / 2
  }

  # -----------------------------
  # 3) g(r) and gradient wrt (rho_xy, rho_xz, rho_yz)
  # -----------------------------
  g_fun_local <- function(r, rho_xy, rho_xz, rho_yz) {
    rho_xy <- clamp(rho_xy)
    r      <- clamp(r)
    rho_xz * r - (rho_xy * rho_xz - rho_yz) * sqrt((1 - r^2) / (1 - rho_xy^2))
  }

  g_grad_local <- function(r, rho_xy, rho_xz, rho_yz) {
    rho_xy <- clamp(rho_xy)
    r      <- clamp(r)
    S <- sqrt((1 - r^2) / (1 - rho_xy^2))
    c(
      # d g / d rho_xy
      -rho_xz * S + (rho_xy * rho_xz - rho_yz) * S * rho_xy / (1 - rho_xy^2),
      # d g / d rho_xz
      r - rho_xy * S,
      # d g / d rho_yz
      S
    )
  }

  # -----------------------------
  # 4) Grid search for min/max over r_xu in D
  # -----------------------------
  D_grid <- seq(rxu_range[1], rxu_range[2], length.out = grid_length)
  D_grid <- clamp(D_grid)

  g_vals <- vapply(D_grid, g_fun_local, numeric(1),
                   rho_xy = rho_xy, rho_xz = rho_xz, rho_yz = rho_yz)

  idx_l <- which.min(g_vals)
  idx_u <- which.max(g_vals)

  r_l <- D_grid[idx_l]
  r_u <- D_grid[idx_u]

  plug_in <- c(g_vals[idx_l], g_vals[idx_u])   # [min, max]

  # -----------------------------
  # 5) Delta-method SE at argmin/argmax + Bonferroni/union-bound critical value
  # -----------------------------
  grad_l <- g_grad_local(r_l, rho_xy, rho_xz, rho_yz)
  grad_u <- g_grad_local(r_u, rho_xy, rho_xz, rho_yz)

  var_l <- as.numeric(t(grad_l) %*% Sigma_rho %*% grad_l)
  var_u <- as.numeric(t(grad_u) %*% Sigma_rho %*% grad_u)

  var_l <- max(var_l, 0)
  var_u <- max(var_u, 0)

  K <- length(D_grid)
  # union bound across r-grid and both endpoints, two-sided:
  crit <- qnorm(1 - alpha / (4 * K))

  CI_lower <- plug_in[1] - crit * sqrt(var_l)
  CI_upper <- plug_in[2] + crit * sqrt(var_u)

  list(
    CI        = c(CI_lower, CI_upper),
    plug_in   = plug_in,
    r_star    = c(r_l, r_u),
    rho_hat   = deltahat,
    Sigma_rho = Sigma_rho,
    crit      = crit
  )
}
