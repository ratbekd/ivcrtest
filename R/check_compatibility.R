#==========================================================
# The main function
#==============================================================
#' Computes the CIs and probabilities
#' @param df Data frame
#' @param n Numeric scalar for size of the variables.
#' @param i Numeric scalar for seed.
#' @param k Numeric scalar for direction of endogenity.
#' @param alpha Numeric scalar for the level of H0 test.
#' @param rxu_range Numeric interval for r_xu values.
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
check_compatibility <- function(df, i = 1, n, k = 1,
                                alpha     = 0.05,
                                rxu_range = rxu_range) {

    # ---- validate inputs ----
  if (!all(c("x", "y", "z") %in% names(df))) {
    stop("DataFrame must contain columns 'x', 'y', 'z'")
  }

  df <- df[complete.cases(df[, c("x", "y", "z")]), ]
  n  <- nrow(df)

  if (n < 30) {
    warning(paste("Small sample size:", n, "observations"))
    return(NULL)
  }

  # ---- sample correlations ----
  rho_xz <- cor(df$x, df$z, use = "complete.obs", method = "pearson")
  rho_yz <- cor(df$y, df$z, use = "complete.obs", method = "pearson")
  rho_xy <- cor(df$y, df$x, use = "complete.obs", method = "pearson")

  # ---- covariance of correlation vector (bootstrap) ----
  deltaSigma <- estimate_cov_corr_boot(df$x, df$y, df$z, B = 800, seed = 1000 + i)
  deltahat   <- c(rho_xy, rho_xz, rho_yz)

  # ---- grid and g function ----
  r_grid <- seq(rxu_range[1], rxu_range[2], length.out = 50)
  g      <- function(delta) g_xu_safe(r_grid, delta)

  # ---- MCUB tuning parameters ----
  B      <- 500
  Blarge <- B * 10
  eta    <- 0.001
  alphac <- 0.8 * alpha
  tol    <- 1e-3
  tol_r  <- 1e-3

  # ---- Jacobian ----
  A  <- t(sapply(r_grid, function(r_xu) {
    as.numeric(local_compute_gradient_safe(r_xu, deltahat[1], deltahat[2], deltahat[3]))
  }))
  Al <- A
  Au <- A

  # ---- g_fun for IM/Stoye (scalar interface) ----
  g_fun_IM <- function(rho_vec, r_xu) {
    rho_xy <- rho_vec["rho_xy"]
    rho_xz <- rho_vec["rho_xz"]
    rho_yz <- rho_vec["rho_yz"]
    rho_xz * r_xu - (rho_xy * rho_xz - rho_yz) *
      sqrt((1 - r_xu^2) / (1 - rho_xy^2))
  }

  # ---- initialise IM outputs (fallback values) ----
  CI_IM          <- c(NA_real_, NA_real_)
  reject_H0_IM   <- NA

  tryCatch({

    cat(sprintf("  Instrument %d: Computing simple CI...\n", i))

    # --- simple union-bound CI ---
    res_simple      <- ci_simple_union(df$x, df$y, df$z,
                                       rxu_range  = rxu_range,
                                       alpha      = alpha,
                                       cov_method = "bootstrap")
    CI_s            <- res_simple$CI
    plug_in         <- res_simple$plug_in
    contains_zero_s <- (CI_s[1] <= 0 & CI_s[2] >= 0)

    cat(sprintf("  Instrument %d: Simple CI computed successfully\n", i))

    # --- MCUB CI (with fallback to simple) ---
    CI_b            <- CI_s
    contains_zero_b <- contains_zero_s
    p_zero          <- ifelse(contains_zero_b, 1, 0.01)

    tryCatch({
      cat(sprintf("  Instrument %d: Computing MCUB CI...\n", i))

      res_bei <- CIhybrid(deltahat, deltaSigma, Al, Au,
                          alpha, alphac, eta, B, Blarge,
                          tol, tol_r, index = NULL, g = g)

      CI_b            <- res_bei$CI_h
      contains_zero_b <- (CI_b[1] <= 0 & CI_b[2] >= 0)

      p_zero <- pvalue_mcub_zero_fast(deltahat, deltaSigma, Al, Au, g,
                                      eta = eta, tol = tol, tol_r = tol_r,
                                      B_fast = 300, Blarge_fast = 3000)

      cat(sprintf("  Instrument %d: MCUB computed successfully\n", i))

    }, error = function(e) {
      cat(sprintf("  Instrument %d: MCUB failed, using simple CI: %s\n",
                  i, e$message))
    })

    # --- IM/Stoye CI ---
    tryCatch({
      cat(sprintf("  Instrument %d: Computing IM/Stoye CI...\n", i))

      res_IM <- im_stoye_inference(
        data      = df,
        grid_r_xu = r_grid,
        g_fun     = g_fun_IM,
        x_col     = "x",
        y_col     = "y",
        z_col     = "z",
        se_method = "bootstrap",
        B         = 999L,
        alpha     = alpha,
        seed      = 1000L + i
      )

      CI_IM        <- res_IM$CI_IM          # named c(lower=..., upper=...)
      reject_H0_IM <- res_IM$reject_H0_0_in_C

      cat(sprintf("  Instrument %d: IM/Stoye CI computed successfully\n", i))

    }, error = function(e) {
      cat(sprintf("  Instrument %d: IM/Stoye failed: %s\n", i, e$message))
    })

    # ---- return results ----
    data.frame(
      Z          = paste0("Z", i),
      plug_in    = sprintf("[%.3f, %.3f]", plug_in[1], plug_in[2]),
      CI_Bei     = sprintf("[%.3f, %.3f]", CI_b[1], CI_b[2]),
      CI_simple  = sprintf("[%.3f, %.3f]", max(-1, CI_s[1]), min(1, CI_s[2])),
      CI_IM      = sprintf("[%.3f, %.3f]",
                           ifelse(is.na(CI_IM[1]), NA_real_, CI_IM["lower"]),
                           ifelse(is.na(CI_IM[2]), NA_real_, CI_IM["upper"])),
      Zero_in_CI_MCUB = ifelse(contains_zero_b,    "\u2713", "\u00d7"),
      Zero_in_CI_IM   = ifelse(isTRUE(!reject_H0_IM), "\u2713", "\u00d7"),
      p_zero     = round(p_zero, 3),
      n          = n,
      stringsAsFactors = FALSE
    )

  }, error = function(e) {
    cat(sprintf("  Instrument %d: Complete failure: %s\n", i, e$message))
    return(NULL)
  })
}
