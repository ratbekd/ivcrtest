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
check_compatibility <- function(data, n, k=-1,
                                alpha = 0.05,
                                seed = 123,
                                rxu_range = c(0, 0.8),
                                bias_mc = TRUE,
                                mc_B = 500) {
  library(dplyr)
  library(ivreg)
  library(lmtest)
  library(sandwich)
  library(AER)
  if (!is.data.frame(data)) stop("Argument 'dataset' must be a data frame.")
  if (!all(c("x", "y", "z") %in% names(data))) {
    stop("Dataset must contain columns named x, y, and z.")
  }

  set.seed(seed)

  rho_xz<-cor(data$x,data$z,use = "complete.obs",method=c("pearson"))
  rho_yz<-cor(data$y,data$z,use = "complete.obs",method=c("pearson"))
  rho_xy<-cor(data$y,data$x,use = "complete.obs",method=c("pearson"))
  # --- Asymptotic covariance matrix under bivariate normality ---
  var_rho_xy <- (1 - rho_xy^2)^2 / n
  var_rho_xz <- (1 - rho_xz^2)^2 / n
  var_rho_yz <- (1 - rho_yz^2)^2 / n
  deltaSigma <- diag(c(var_rho_xy, var_rho_xz, var_rho_yz))

  # --- Point estimates ---
  deltahat <- c(rho_xy, rho_xz, rho_yz)

  # --- Grid of possible r_xu values ---
  r_grid <- seq(rxu_range[1], rxu_range[2], length.out = 50)
  D_grid = seq(rxu_range[1], rxu_range[2], length.out = 50)
  # --- g(delta) mapping ---
  #g <- function(delta) g_xu_safe(r_grid, matrix(delta, nrow = 1))
  g <- function(delta) g_xu_safe(r_grid, delta)

  # --- Tuning parameters for Union Bound ---
  B      <- 500
  Blarge <- B * 10
  eta    <- 0.001
  alphac <- 0.8 * alpha
  tol    <- 1e-3
  tol_r  <- 1e-3

  # --- Compute A matrix (Jacobian of g wrt delta) ---
  A <- t(sapply(r_grid, function(r_xu) {
    grad <- local_compute_gradient_safe(r_xu, deltahat[1], deltahat[2], deltahat[3])
    as.numeric(grad)
  }))
  Al <- A
  Au <- A

  # Compute simple union CI
  res_simple <- ci_simple_union(data$x, data$y, data$z,
                                rxu_range=rxu_range,
                                alpha = alpha)


  res_bei <-CIhybrid(deltahat, deltaSigma, Al, Au, alpha, alphac,
                     eta, B, Blarge, tol, tol_r, index = NULL, g = g)

  CI_b <- res_bei$CI_h
  # CI_s <- CI_simple
  CI_s <- res_simple$CI
  # -------- g_fun plug-in CI --------

  plug <- res_bei$CI_c
  plug_in <- res_simple$plug_in
  contains_zero_b <- (CI_b[1] <= 0 & CI_b[2] >= 0)
  contains_zero_s <- (CI_s[1] <= 0 & CI_s[2] >= 0)

  # coverage test: does CI cover plug-in interval
  cover_b <- (plug[1] >= CI_b[1] & plug[2] <= CI_b[2])
  cover_s <- (plug[1] >= CI_s[1] & plug[2] <= CI_s[2])
  ### === 9. CI for r_zu = 0 ===
  grad_rzu0 <- g_grad(0, rho_xy, rho_xz, rho_yz)

  # rho_xu function
  r_xy <- cor(data$x, data$y)
  r_xz <- cor(data$x, data$z)
  r_yz <- cor(data$y, data$z)
  r_xz<-cor(data$x,data$z,use = "complete.obs",method=c("pearson"))
  r_yz<-cor(data$y,data$z,use = "complete.obs",method=c("pearson"))
  r_xy<-cor(data$y,data$x,use = "complete.obs",method=c("pearson"))
  rho_xu <- function(r_xz, r_xy, r_yz) {
    mult <- 1 - r_xy^2
    if (mult <= 1e-10) return(NA)
    denom <- (r_xy * r_xz - r_yz)^2
    if (denom <= 1e-10) return(NA)
    num <- pmax(0, r_xz^2)
    res <- sqrt(1 / (1 + mult * num / denom))
    if (!is.finite(res)) return(NA)
    res
  }

  rxu_point <- k*rho_xu(r_xz, r_xy, r_yz)
  data <- data
  n=200
  # Delta method for rxu_point
  delta_method_rho_xu <- function(data, r_xy, r_xz, r_yz, n,
                                  grad_rzu0, rxu_point, alpha = 0.05) {
    if (is.null(data) || is.na(rxu_point)) {
      return(list(
        point_estimate = rxu_point,
        point_estimate_bias_corrected = rxu_point,
        ci_bias_corrected = c(L = NA, U = NA),
        bias_correction = 0
      ))
    }

    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      se_simple <- 1/sqrt(n)
      z <- qnorm(1 - alpha/2)
      return(list(
        point_estimate = rxu_point,
        point_estimate_bias_corrected = rxu_point,
        ci_bias_corrected = c(L = rxu_point - z*se_simple,
                              U = rxu_point + z*se_simple),
        bias_correction = 0
      ))
    }

    theta_hat <- c(r_xy, r_xz, r_yz)

    g_fun_local <- function(theta) {
      theta[2] * rxu_point - (theta[1] * theta[2] - theta[3]) * sqrt((1 - rxu_point^2) / (1 - theta[1]^2))
    }

    grad_g0 <- grad_rzu0

    boot_corrs <- replicate(200, {
      idx <- sample(n, replace = TRUE)
      d <- data[idx, ]
      c(cor(d$x, d$y), cor(d$x, d$z), cor(d$y, d$z))
    })
    Sigma_hat <- cov(t(boot_corrs))

    var_g <- as.numeric(t(grad_g0) %*% Sigma_hat %*% grad_g0)
    se_g <- sqrt(var_g / n)

    hess_g <- numDeriv::hessian(g_fun_local, theta_hat)
    bias <- sum(0.5 * hess_g * Sigma_hat) / n
    rxu_bc <- rxu_point - bias

    z <- qnorm(1 - alpha/2)
    ci <- c(L = rxu_bc - z*se_g, U = rxu_bc + z*se_g)

    list(
      point_estimate = rxu_point,
      point_estimate_bias_corrected = rxu_bc,
      ci_bias_corrected = ci,
      bias_correction = bias
    )
  }

  rxuf <- delta_method_rho_xu(data, r_xy, r_xz, r_yz, n,
                              grad_rzu0, rxu_point, alpha = 0.05)

  rxu_point_corrected <- rxuf$point_estimate_bias_corrected
  ci_rxu <- rxuf$ci_bias_corrected

  # Bias correction for target value
  bias_z <- 0
  range_interval <- NULL
  bias_mc <- TRUE
  mc_B <- 200

  if (bias_mc && !is.null(data)) {
    boot_corrs <- replicate(mc_B, {
      idx <- sample(nrow(data), replace = TRUE)
      d <- data[idx, ]
      c(xz = cor(d$x, d$z, use = "complete.obs"),
        xy = cor(d$x, d$y, use = "complete.obs"),
        yz = cor(d$y, d$z, use = "complete.obs"))
    })
    valid <- complete.cases(t(boot_corrs))
    boot_corrs <- boot_corrs[, valid]

    if (ncol(boot_corrs) > 100) {
      sim_bias_z <- apply(boot_corrs, 2, function(corrs) {
        g_fun(corrs["xz"], corrs["xy"], corrs["yz"], rxu_point_corrected) -
          g_fun(r_xz, r_xy, r_yz, rxu_point_corrected)
      })
      bias_z <- mean(sim_bias_z, na.rm = TRUE)
      range_interval <- quantile(sim_bias_z, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
    }
  }

  bias_corrected_g <- g_fun(rxu_point_corrected, r_xy, r_xz, r_yz) - bias_z
  bias_corrected_g <-  bias_z
  # bias_corrected_g <- g_fun(rxu_point, r_xy, r_xz, r_yz)
  target_v <- bias_corrected_g
  target_interval <- c(bias_corrected_g + range_interval[1],
                       bias_corrected_g + range_interval[2])

  CI_rzu0 <- target_interval
  cat(sprintf("CI for r_zu = 0: [%.4f, %.4f]\n", CI_rzu0[1], CI_rzu0[2]))
  CI_rxu_z0=c(g_fun(ci_rxu[1], r_xy, r_xz, r_yz),g_fun(ci_rxu[2], r_xy, r_xz, r_yz))

  #

  target_v <- if(CI_rzu0[1]*CI_rzu0[2]>0) 0 else target_v
  # P-value for zero in CI
  # Approximate SE from CI width
  se_hat <- (CI_b[2] - CI_b[1]) / (2 * qnorm(0.975))

  # z-stat: how far target_v is outside CI (0 if inside)
  if (target_v < CI_b[1]) {
    z_stat <- (CI_b[1] - target_v) / se_hat
  } else if (target_v > CI_b[2]) {
    z_stat <- (target_v - CI_b[2]) / se_hat
  } else {
    z_stat <- 0
  }

  # Two-sided p-value; if inside CI, set p=1
  p_zero <- if(z_stat == 0) 1 else 2 * (1 - pnorm(abs(z_stat)))

  ### === 10. Overlap and containment tests ===
  # Full function to compute Sigma_g
  compute_Sigma_g <- function(data, D_grid, g_grad, g_fun, method = "bootstrap") {

    # Compute correlations
    # rho_xy <- cor(data$x, data$y)
    # rho_xz <- cor(data$x, data$z)
    # rho_yz <- cor(data$y, data$z)
    rho_xz<-cor(data$x,data$z,use = "complete.obs",method=c("pearson"))
    rho_yz<-cor(data$y,data$z,use = "complete.obs",method=c("pearson"))
    rho_xy<-cor(data$y,data$x,use = "complete.obs",method=c("pearson"))
    n <- 200

    # Find min/max and their achieving r values
    g_vals <- sapply(D_grid, function(r) g_fun(r, rho_xy, rho_xz, rho_yz))
    g_min <- min(g_vals, na.rm = TRUE)
    g_max <- max(g_vals, na.rm = TRUE)

    idx_min <- which.min(g_vals)
    idx_max <- which.max(g_vals)
    r_min <- D_grid[idx_min]
    r_max <- D_grid[idx_max]

    if (method == "delta") {
      # Asymptotic variance of correlations
      var_xy <- (1 - rho_xy^2)^2 / n
      var_xz <- (1 - rho_xz^2)^2 / n
      var_yz <- (1 - rho_yz^2)^2 / n
      Sigma_rho <- diag(c(var_xy, var_xz, var_yz))

      # Gradients
      grad_min <- g_grad(r_min, rho_xy, rho_xz, rho_yz)
      grad_max <- g_grad(r_max, rho_xy, rho_xz, rho_yz)

      # Delta method
      Jacobian <- rbind(grad_min, grad_max)
      Sigma_g <- Jacobian %*% Sigma_rho %*% t(Jacobian)

    } else if (method == "bootstrap") {
      # Bootstrap correlations
      B <- 500
      boot_corrs <- replicate(B, {
        idx <- sample(n, replace = TRUE)
        d <- data[idx, ]
        c(cor(d$x, d$y), cor(d$x, d$z), cor(d$y, d$z))
      })
      Sigma_rho <- cov(t(boot_corrs))

      # Gradients at observed r values
      grad_min <- g_grad(r_min, rho_xy, rho_xz, rho_yz)
      grad_max <- g_grad(r_max, rho_xy, rho_xz, rho_yz)

      # Delta method with bootstrap Sigma_rho
      Jacobian <- rbind(grad_min, grad_max)
      Sigma_g <- Jacobian %*% Sigma_rho %*% t(Jacobian)

    } else if (method == "direct_bootstrap") {
      # Bootstrap (g_min, g_max) directly
      B <- 500
      boot_bounds <- replicate(B, {
        idx <- sample(n, replace = TRUE)
        d <- data[idx, ]

        rho_xy_b <- cor(d$x, d$y)
        rho_xz_b <- cor(d$x, d$z)
        rho_yz_b <- cor(d$y, d$z)

        g_vals_b <- sapply(D_grid, function(r) g_fun(r, rho_xy_b, rho_xz_b, rho_yz_b))

        c(min(g_vals_b, na.rm = TRUE), max(g_vals_b, na.rm = TRUE))
      })

      Sigma_g <- cov(t(boot_bounds))
    }

    return(list(
      Sigma_g = Sigma_g,
      r_min = r_min,
      r_max = r_max,
      g_min = g_min,
      g_max = g_max
    ))
  }

  # Usage in your loop:
  result <- compute_Sigma_g(data, D_grid, g_grad, g_fun, method = "delta")
  Sigma_g <- result$Sigma_g
  # Function to compute coverage probability
  #compute_plugin_coverage_prob <- function(CI, plug_in=plug_in_set, Sigma_g, method = "simulation") {
    compute_plugin_coverage_prob <- function(CI, plug_in, Sigma_g, method = "simulation"){

    if (method == "simulation") {
      B_sim <- 5000
      g_sim <- MASS::mvrnorm(B_sim, mu = plug_in, Sigma = Sigma_g)
      contains <- (CI[1] <= g_sim[,1]) & (CI[2] >= g_sim[,2])
      return(mean(contains))

    } else if (method == "bootstrap") {
      # Requires data - see Method 1 above
      stop("Bootstrap method requires data input")

    } else if (method == "normal") {
      # Check if CI currently contains plug-in
      if (CI[1] <= plug_in[1] && CI[2] >= plug_in[2]) {
        # Compute probability it remains contained
        lower_gap <- plug_in[1] - CI[1]
        upper_gap <- CI[2] - plug_in[2]

        se_lower <- sqrt(Sigma_g[1,1])
        se_upper <- sqrt(Sigma_g[2,2])

        z_lower <- lower_gap / se_lower
        z_upper <- upper_gap / se_upper

        return(pnorm(z_lower) * pnorm(z_upper))
      } else {
        # Compute probability of containment
        lower_gap <- abs(plug_in[1] - CI[1])
        upper_gap <- abs(CI[2] - plug_in[2])

        se_lower <- sqrt(Sigma_g[1,1])
        se_upper <- sqrt(Sigma_g[2,2])

        z_lower <- lower_gap / se_lower
        z_upper <- upper_gap / se_upper

        if (CI[1] > plug_in[1] && CI[2] < plug_in[2]) {
          # CI is inside plug-in set
          return(1 - (pnorm(-z_lower) + pnorm(-z_upper)))
        } else {
          return((1 - pnorm(z_lower)) * (1 - pnorm(z_upper)))
        }
      }
    }
  }
  overlap <- !(CI_rzu0[2] < CI_b[1] | CI_rzu0[1] > CI_b[2])
  contain <- (CI_rzu0[1] >= CI_b[1]  & CI_rzu0[2] <= CI_b[2] )
  contains_zero_TI <- (CI_rzu0[1] <= 0 & CI_rzu0[2] >= 0)

  # #Pvalue of coverage of target interval by CI
 #

  p_contain_s <- compute_plugin_coverage_prob(CI=CI_s, plug_in=CI_rzu0, Sigma_g, method = "simulation")

  #############################Experiment with coverage
  # Compute Sigma_g for coverage probability calculation
  # Need to compute covariance matrix of (g_min, g_max)
  n <- round(nrow(data))
  D_grid <- seq(0.00, 0.80, length.out = 50)

  # Find r values that achieve min/max
  g_vals <- sapply(D_grid, function(r) g_fun(r, rho_xy, rho_xz, rho_yz))
  idx_min <- which.min(g_vals)
  idx_max <- which.max(g_vals)
  r_min <- D_grid[idx_min]
  r_max <- D_grid[idx_max]

  # Compute gradients at min/max achieving r values
  grad_min <- g_grad(r_min, rho_xy, rho_xz, rho_yz)
  grad_max <- g_grad(r_max, rho_xy, rho_xz, rho_yz)

  # Variance estimates
  var_xy <- (1 - rho_xy^2)^2 / n
  var_xz <- (1 - rho_xz^2)^2 / n
  var_yz <- (1 - rho_yz^2)^2 / n
  Sigma_rho <- diag(c(var_xy, var_xz, var_yz))

  # Covariance of (g_min, g_max) using delta method
  Jacobian <- rbind(grad_min, grad_max)
  Sigma_g <- Jacobian %*% Sigma_rho %*% t(Jacobian)

  # Now compute coverage probability
  p_coverage <- compute_plugin_coverage_prob(CI_b, plug_in, Sigma_g, method = "normal")

  # Check if zero is in CI
  contains_zero_b <- (CI_b[1] <= 0 & CI_b[2] >= 0)

  # Check if CI covers plug-in interval
  cover_b <- (plug_in[1] >= CI_b[1] & plug_in[2] <= CI_b[2])

  # p value for coverage of plug-in by CI
  p_coverage_s <- compute_plugin_coverage_prob(CI_s, plug_in, Sigma_g, method = "normal")

  #Pvalue of coverage of target interval by CI
  p_contain <- compute_plugin_coverage_prob(CI=CI_b, plug_in=CI_rzu0, Sigma_g, method = "simulation")
  #############################################################

  # Variance-covariance matrix of correlations (asymptotic)
  var_rxy <- (1 - r_xy^2)^2 / n
  var_rxz <- (1 - r_xz^2)^2 / n
  var_ryz <- (1 - r_yz^2)^2 / n

  # More realistic correlation structure
  if (!is.null(data)) {
    # Use bootstrap to estimate correlation covariance
    boot_corrs <- replicate(500, {
      idx <- sample(nrow(data), replace = TRUE)
      d <- data[idx, ]
      c(cor(d$x, d$y), cor(d$x, d$z), cor(d$y, d$z))
    })
    vcov_rhos <- cov(t(boot_corrs))
  } else {
    # Simplified assumption of independence
    vcov_rhos <- diag(c(var_rxy, var_rxz, var_ryz))
  }
  rxu_range=c(min(D_grid), max(D_grid))
  # Delta method standard errors at optimal points
  grad_min <- g_grad(r_xy, r_xz, r_yz, rxu_range[1])
  grad_max <- g_grad(r_xy, r_xz, r_yz, rxu_range[2])

  if (any(is.na(grad_min)) || any(is.na(grad_max))) {
    se_min <- sqrt(var(g_vals) / n)  # fallback
    se_max <- sqrt(var(g_vals) / n)  # fallback
  } else {
    se_min <- sqrt(as.numeric(t(grad_min) %*% vcov_rhos %*% grad_min))
    se_max <- sqrt(as.numeric(t(grad_max) %*% vcov_rhos %*% grad_max))
  }
  ################ function overlap pval calculation
  overlap_test_interval <- function(CI, target_interval, se_min, se_max) {
    # No overlap case
    if (target_interval[2] <= CI[1]) {
      z_stat <- (CI[1] - target_interval[2]) / se_min
      p_val <- 1 - pnorm(z_stat)
      return(list(overlap = FALSE, p_value_overlap = p_val))
    } else if (target_interval[1] >= CI[2]) {
      z_stat <- (target_interval[1] - CI[2]) / se_max
      p_val <- 1 - pnorm(z_stat)
      return(list(overlap = FALSE, p_value_overlap = p_val))
    } else {
      # Overlap exists → high p-value
      return(list(overlap = TRUE, p_value_overlap = 0.95))
    }
  }
  res_overlap <- overlap_test_interval(CI=CI_b, target_interval, se_min, se_max)

  p_overlap = res_overlap$p_value_overlap#,


  data.frame(
    Z = paste0("IV"),
    plug_in = sprintf("[%.2f,%.2f]", plug_in[1],plug_in[2]),
    CI_Bei = sprintf("[%.2f,%.2f]", CI_b[1], CI_b[2]),
    Plug_covered_CI_Bei = ifelse(cover_b, "✓", "×"),
    p_coverage = round(p_coverage, 3),
    CI_simple = sprintf("[%.2f,%.2f]", max(-1,CI_s[1]), min(1,CI_s[2])),
    Plug_covered_CI_s = ifelse(cover_s, "✓", "×"),
    p_coverage_s=round(p_coverage_s, 3),
    Zero_in_CI_Bei = ifelse(contains_zero_b, "✓", "×"),
    p_zero= round(p_zero,3),
    Target_interval=sprintf("[%.2f,%.2f]", CI_rzu0[1], CI_rzu0[2]),
    Target_interval_overlap=ifelse(overlap, "✓", "×"),
    Prob_Target_inter_overlap=round(p_overlap,3),
    CI_Bei_contains_TI=ifelse(contain, "✓", "×"),
    Prob_Target_inter_cont=round(p_contain,3),
    Prob_Target_inter_cont_s=round(p_contain_s,3)

  )
}

