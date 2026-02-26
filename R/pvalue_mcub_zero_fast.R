#' Fast p-value computation for MCUB zero test
#' @param deltahat Numeric vector of sample correlations.
#' @param deltaSigma Covariance matrix of deltahat.
#' @param Al Lower Jacobian matrix.
#' @param Au Upper Jacobian matrix.
#' @param g Function mapping delta to lambda.
#' @param eta Numeric tuning parameter.
#' @param B_fast Integer number of fast bootstrap draws.
#' @param Blarge_fast Integer number of large bootstrap draws.
#' @param tol Numeric tolerance.
#' @param tol_r Numeric tolerance for correlations.
#' @param alpha_grid Numeric vector of alpha values for grid search.
#' @param refine_steps Integer number of bisection refinement steps.
#' @param seed Integer seed for reproducibility.
#' @return Numeric p-value.
#' @export
pvalue_mcub_zero_fast <- function(deltahat, deltaSigma, Al, Au, g,
                                  eta = 0.001,
                                  B_fast = 600, Blarge_fast = 6000,
                                  tol = 1e-3, tol_r = 1e-3,
                                  alpha_grid = c(1e-4, 5e-4, 1e-3, 2e-3, 5e-3,
                                                 0.01, 0.02, 0.05, 0.10, 0.20, 0.40, 0.70),
                                  refine_steps = 6,
                                  seed = 123) {

  inside0 <- function(alpha) {
    set.seed(seed + as.integer(alpha * 1e6))  # stabilize across calls
    alpha_C <- alpha / 2
    alphac  <- 0.8 * alpha_C

    res <- CIhybrid(deltahat, deltaSigma, Al, Au,
                    alpha = alpha_C, alphac = alphac,
                    eta = eta, B = B_fast, Blarge = Blarge_fast,
                    tol = tol, tol_r = tol_r, index = NULL, g = g)
    CI <- res$CI_h
    (CI[1] <= 0 && 0 <= CI[2])
  }

  # 1) Grid search
  ins <- vapply(alpha_grid, inside0, logical(1))

  # If 0 excluded even at tiny alpha -> p ~ 0
  if (!ins[1]) return(alpha_grid[1])

  # If 0 included even at huge alpha -> p ~ 1
  if (ins[length(ins)]) return(1.0)

  # Find first alpha where 0 becomes excluded
  j <- which(!ins)[1]
  lo <- alpha_grid[j - 1]  # included
  hi <- alpha_grid[j]      # excluded

  # 2) Short refinement
  for (k in seq_len(refine_steps)) {
    mid <- 0.5 * (lo + hi)
    if (inside0(mid)) lo <- mid else hi <- mid
  }

  hi
}
