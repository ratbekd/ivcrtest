#==========================================================
# The main function
#==============================================================
#' Computes the CIs and probabilities
#' @param data Data frame
#' @param n Numeric scalar for size of the variables.
#' @param k Numeric scalar for direction of endogenity.
#' @param alpha Numeric scalar for the level of H0 test.
#' @param rxu_range Numeric interval for r_xu values.
#' @param seed Numeric scalar  for reproducibility.
#' @param bias_mc Logial value.
#' @param mc_B Numeric scalar for MC correction simulations.
#' @return CIs and probabilities.
#' @export
#' @importFrom stats  predict
#' @importFrom stats  complete.cases
#' @importFrom stats  fitted
#' @importFrom stats pf
#' @importFrom stats pchisq
#' @importFrom stats var
#' @importFrom stats as.formula
#' @importFrom stats lm
#' @importFrom stats cor
#' @importFrom stats cov
#' @importFrom stats residuals
#' @importFrom stats df.residual
#' @importFrom stats qt
#' @importFrom stats sd
#' @importFrom stats rnorm
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats predict
#' @importFrom stats df
#' @importFrom stats optim
#' @importFrom stats resid
#' @importFrom stats quantile
check_compatibility <-  function(df, i=1, n, k = 1,
                                 alpha = 0.05,
                                 rxu_range = rxu_range) {
  library(dplyr)
  library(ivreg)
  library(lmtest)
  library(sandwich)
  library(AER)
  # Ensure column names are correct
  if (!all(c("x", "y", "z") %in% names(df))) {
    stop("DataFrame must contain columns 'x', 'y', 'z'")
  }

  # Remove any rows with missing values
  df <- df[complete.cases(df[, c("x", "y", "z")]), ]
  n <- nrow(df)

  if (n < 30) {
    warning(paste("Small sample size:", n, "observations"))
    return(NULL)
  }

  # Compute correlations
  rho_xz <- cor(df$x, df$z, use = "complete.obs", method = "pearson")
  rho_yz <- cor(df$y, df$z, use = "complete.obs", method = "pearson")
  rho_xy <- cor(df$y, df$x, use = "complete.obs", method = "pearson")

  # Asymptotic covariance
  var_rho_xy <- (1 - rho_xy^2)^2 / n
  var_rho_xz <- (1 - rho_xz^2)^2 / n
  var_rho_yz <- (1 - rho_yz^2)^2 / n
  deltaSigma <- diag(c(var_rho_xy, var_rho_xz, var_rho_yz))

  deltahat <- c(rho_xy, rho_xz, rho_yz)

  # Grid for r_xu
  r_grid <- seq(rxu_range[1], rxu_range[2], length.out = 50)

  # Define g function
  g <- function(delta) g_xu_safe(r_grid, delta)

  # MCUB parameters
  B <- 500  # Reduced for speed
  Blarge <- B * 10
  eta <- 0.001
  alphac <- 0.8 * alpha
  tol <- 1e-3
  tol_r <- 1e-3

  # Compute Jacobian
  A <- t(sapply(r_grid, function(r_xu) {
    grad <- local_compute_gradient_safe(r_xu, deltahat[1], deltahat[2], deltahat[3])
    as.numeric(grad)
  }))
  Al <- A
  Au <- A

  # Compute simple CI first (faster)
  tryCatch({
    cat(sprintf("  Instrument %d: Computing simple CI...\n", i))

    # Pass vectors, not dataframe
    res_simple <- ci_simple_union(df$x, df$y, df$z,
                                  rxu_range = rxu_range,
                                  alpha = alpha)

    CI_s <- res_simple$CI
    plug_in <- res_simple$plug_in
    contains_zero_s <- (CI_s[1] <= 0 & CI_s[2] >= 0)

    cat(sprintf("  Instrument %d: Simple CI computed successfully\n", i))

    # Try MCUB (this may fail, so we have fallback)
    CI_b <- CI_s  # Default to simple
    contains_zero_b <- contains_zero_s
    p_zero <- ifelse(contains_zero_b, 1, 0.01)

    tryCatch({
      cat(sprintf("  Instrument %d: Computing MCUB CI...\n", i))

      res_bei <- CIhybrid(deltahat, deltaSigma, Al, Au, alpha, alphac,
                          eta, B, Blarge, tol, tol_r, index = NULL, g = g)

      CI_b <- res_bei$CI_h
      contains_zero_b <- (CI_b[1] <= 0 & CI_b[2] >= 0)

      # P-value computation
      p_zero <- pvalue_mcub_zero_fast(deltahat, deltaSigma, Al, Au, g,
                                      eta = eta, tol = tol, tol_r = tol_r,
                                      B_fast = 300, Blarge_fast = 3000)

      cat(sprintf("  Instrument %d: MCUB computed successfully\n", i))

    }, error = function(e) {
      cat(sprintf("  Instrument %d: MCUB failed, using simple CI: %s\n", i, e$message))
    })

    # Return results
    data.frame(
      Z = paste0("Z", i),
      # r_xz = round(rho_xz, 3),
      # r_zu = if ("u" %in% names(df)) round(cor(df$z, df$u), 3) else NA,
      # beta_IV = round(cov(df$z, df$y, use = "complete.obs") /
      #                   cov(df$z, df$x, use = "complete.obs"), 3),
      plug_in = sprintf("[%.3f, %.3f]", plug_in[1], plug_in[2]),
      CI_Bei = sprintf("[%.3f, %.3f]", CI_b[1], CI_b[2]),
      CI_simple = sprintf("[%.3f, %.3f]", max(-1, CI_s[1]), min(1, CI_s[2])),
      Zero_in_CI = ifelse(contains_zero_b, "✓", "×"),
      p_zero = round(p_zero, 3),
      n = n,
      stringsAsFactors = FALSE
    )

  }, error = function(e) {
    cat(sprintf("  Instrument %d: Complete failure: %s\n", i, e$message))
    return(NULL)
  })
}


