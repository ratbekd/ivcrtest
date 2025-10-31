#'  Returns final confidence intervals
#' @param deltahat row or column vectors/matrices.
#' @param Al constraint matrices
#' @param Au constraint matrices
#' @param sigma_l standard deviations of lambda_hat_l / u
#' @param sigma_u standard deviations of lambda_hat_l / u
#' @param deltaSigma (B x k) matrix
#' @param alpha,  critical levels
#' @param alphac critical levels
#' @param eta critical area level
#' @param index index number
#' @param B,  number of repetitions for bootstrap
#' @param Blarge number of repetitions for bootstrap
#' @param tol tolerance for correlation threshold
#' @param tol_r tolerance for correlation threshold
#' @param g mapping function returning lambda estimates
#' @return Return final confidence intervals
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
#'
CIhybrid <- function(deltahat, deltaSigma, Al, Au,
                     alpha, alphac, eta,
                     B, Blarge,
                     tol, tol_r, index, g) {

  set.seed(0)

  # ----- 1. Random draws and basic quantities -----
  library(MASS)
  library(stats)

  k_dim <- length(deltahat)
  deltastar_demean_large <- MASS::mvrnorm(Blarge, mu = rep(0, k_dim), Sigma = deltaSigma)
  deviation <- apply(abs(deltastar_demean_large / sqrt(diag(deltaSigma))), 1, max)
  c_bd <- quantile(deviation, 1 - eta)

  sigma_delta <- sqrt(diag(deltaSigma))

  # ----- 2. Construct lambdaSigma and related components -----
  Lambda <- rbind(Al, Au)
  lambdaSigma <- Lambda %*% deltaSigma %*% t(Lambda)
  lambdasigma <- sqrt(diag(lambdaSigma))

  kk <- nrow(Al)
  sigma_l <- lambdasigma[1:kk]
  sigma_u <- lambdasigma[(kk + 1):(2 * kk)]

  corr_all <- diag(1 / lambdasigma) %*% lambdaSigma %*% diag(1 / lambdasigma)
  corr_l <- corr_all[1:kk, 1:kk]
  corr_u <- corr_all[(kk + 1):(2 * kk), (kk + 1):(2 * kk)]
  corr_m <- corr_all[1:kk, (kk + 1):(2 * kk)]

  # ----- 3. Feasible bounds -----
  lb <- deltahat - sigma_delta * c_bd
  ub <- deltahat + sigma_delta * c_bd
  delta1 <- deltahat

  lb[index] <- 0
  ub[index] <- 0
  delta1[index] <- 0

  cl <- 0
  cu <- qnorm(1 - alpha / 2)
  c <- (cl + cu) / 2

  delta_fea <- list()
  c_fea <- numeric()

  # ----- 4. Objective wrapper -----
  obj_large <- function(delta, c_check) {
    (alpha - CIproj_p(c_check, delta, alphac, c_bd,
                      deltastar_demean_large,
                      Al, Au, sigma_l, sigma_u,
                      deltaSigma, corr_m, corr_l, corr_u,
                      eta, tol_r, g)) * 100
  }

  # ----- 5. Iterative bisection loop -----
  K1 <- 1
  while ((cu - cl) > tol) {
    set.seed(K1)
    K1 <- K1 + 1

    deltastar_demean <- MASS::mvrnorm(B, mu = rep(0, k_dim), Sigma = deltaSigma)

    obj <- function(delta) {
      (alpha - CIproj_p(c, delta, alphac, c_bd, deltastar_demean,
                        Al, Au, sigma_l, sigma_u, deltaSigma,
                        corr_m, corr_l, corr_u, eta, tol_r, g)) * 100
    }

    p1 <- obj_large(delta1, c)

    if (!is.finite(p1)) {
      cat("Warning: p1 is non-finite. delta1 =", round(delta1, 3), " c =", round(c, 3), "\n")
    }

    if (p1 >= 0) {
      # Optimizer equivalent of fmincon + GlobalSearch
      f_optim <- tryCatch({
        optim(par = delta1, fn = obj, method = "L-BFGS-B",
              lower = lb, upper = ub, control = list(maxit = 1000))
      }, error = function(e) NULL)

      if (!is.null(f_optim)) {
        delta1 <- f_optim$par
      }
      p1 <- obj_large(delta1, c)

      if (p1 >= 0) {
        cu <- c
      } else {
        cl <- c
        delta_fea[[length(delta_fea) + 1]] <- delta1
        c_fea <- c(c, c_fea)
      }
    } else {
      cl <- c
      delta_fea[[length(delta_fea) + 1]] <- delta1
      c_fea <- c(c, c_fea)
    }

    c <- (cl + cu) / 2
  }

  # ----- 6. Final refinement -----
  cl <- c
  cu <- qnorm(1 - alpha / 2)

  while ((cu - cl) > tol) {
    c <- (cl + cu) / 2
    p <- obj_large(delta1, c)
    if (p >= 0) {
      cu <- c
    } else {
      cl <- c
    }
  }

  # ----- 7. Compute confidence intervals -----
  lambdahat_l <- g(deltahat)
  lambdahat_u <- lambdahat_l

  c_LF <- qnorm(1 - eta)
  CI_p <- c(min(lambdahat_l - c * sigma_l), max(lambdahat_u + c * sigma_u))
  CI_c <- CIcon(deltahat, deltaSigma, Al, Au, c_LF, alphac, tol, tol_r, g)

  CI_h <- c(min(CI_p[1], CI_c[1]), max(CI_p[2], CI_c[2]))

  list(CI_h = CI_h,
       CI_c = CI_c,
       CI_p = CI_p,
       delta1 = delta1,
       c = c)
}
