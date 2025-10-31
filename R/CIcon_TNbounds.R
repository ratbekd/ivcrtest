#' Computes list(th_1, th_2): lower / upper bounds vectors
#' @param theta scalar
#' @param deltahat (B x k) matrix
#' @param Al constraint matrice
#' @param Au constraint matrices
#' @param sigma_l standard deviations of lambda_hat_l / u
#' @param sigma_u standard deviations of lambda_hat_l / u
#' @param corr_m correlation blocks
#' @param corr_l correlation blocks
#' @param corr_u correlation blocks
#' @param cLFl critical values
#' @param cLFu critical values
#' @param tol_r tolerance for correlation threshold
#' @param g mapping function returning lambda estimates
#' @return list(th_1, th_2): lower / upper bounds vectors
#' @export
CIcon_TNbounds <- function(theta, deltahat,
                           Al, Au,
                           sigma_l, sigma_u,
                           corr_m, corr_l, corr_u,
                           cLFl, cLFu, tol_r, g) {
  # -----------------------------------------------
  #
  # -----------------------------------------------
  # Inputs:
  #  - theta: scalar
  #  - deltahat: (B x k) matrix
  #  - Al, Au: constraint matrices
  #  - sigma_l, sigma_u: standard deviations of lambda_hat_l / u
  #  - corr_m, corr_l, corr_u: correlation blocks
  #  - cLFl, cLFu: critical values
  #  - tol_r: tolerance for correlation threshold
  #  - g: mapping function returning lambda estimates
  #
  # Outputs:
  #  - list(th_1, th_2): lower / upper bounds vectors

  B <- nrow(deltahat)
  k_l <- nrow(Al)
  k_u <- nrow(Au)

  # ----- 1. Lambda estimates -----
  lambdahat_l <- g(deltahat)
  lambdahat_u <- lambdahat_l  # same in original MATLAB

  # ----- 2. Compute T statistics -----
  Tlb <- sweep(lambdahat_l - theta, 2, sigma_l, "/")
  Tub <- sweep(theta - lambdahat_u, 2, sigma_u, "/")

  Tl <- apply(Tlb, 1, min)
  bl <- apply(Tlb, 1, which.min)
  Tu <- apply(Tub, 1, min)
  bu <- apply(Tub, 1, which.min)

  # Extract relevant correlation rows/cols
  #corr_m_blb <- corr_m[cbind(bl, 1:ncol(corr_m))]  # this will be reshaped
  corr_m_blb <- corr_m[bl, , drop = FALSE]
  corr_m_bub <- t(corr_m[, bu, drop = FALSE])

  # ----- 3. Case: Tl >= Tu (lower side) -----
  # TS part
  tTS_l <- matrix(1e10, nrow = B, ncol = k_u)
  tTS_ll <- (1 + corr_m_blb)^(-1) * (matrix(Tub, nrow = B, ncol = k_u) + corr_m_blb * Tl)
  tTS_l[(1 + corr_m_blb) > tol_r] <- tTS_ll[(1 + corr_m_blb) > tol_r]

  th_1_1 <- apply(tTS_l, 1, min) * (Tl >= Tu)
  th_1_1[th_1_1 == 1e10] <- -1e10
  th_1_1 <- th_1_1 * (Tl >= Tu)

  # LF part
  tLFl <- rep(cLFl, B)

  # B part
  # taux_1 <- (1 - corr_l[cbind(bl, 1:ncol(corr_l))])^(-1) *
  #   (Tlb[cbind(1:B, bl)] - corr_l[cbind(bl, 1:ncol(corr_l))] * Tl)

  corr_l_rows <- corr_l[cbind(bl, bl)]
  taux_1 <- (1 - corr_l_rows)^(-1) * (Tlb[cbind(1:B, bl)] - corr_l_rows * Tl)

  tB_l2 <- matrix(1e10, nrow = B, ncol = k_l)
  mask1 <- (1 > corr_l[bl, , drop = FALSE] + tol_r)
  tB_l2[mask1] <- taux_1[mask1]
  th_2_1 <- apply(cbind(tLFl, tB_l2), 1, min) * (Tl >= Tu)

  # ----- 4. Case: Tl < Tu (upper side) -----
  # TS part
  tTS_u <- matrix(1e10, nrow = B, ncol = k_l)
  tTS_u1 <- (1 + corr_m_bub)^(-1) * (matrix(Tlb, nrow = B, ncol = k_l) + corr_m_bub * Tu)
  tTS_u[(1 + corr_m_bub) > tol_r] <- tTS_u1[(1 + corr_m_bub) > tol_r]

  th_1_2 <- apply(tTS_u, 1, min) * (Tl < Tu)
  th_1_2[th_1_2 == 1e10] <- -1e10
  th_1_2 <- th_1_2 * (Tl < Tu)

  # LF part
  tLFu <- rep(cLFu, B)

  # B part
  #taux_2 <- (1 - corr_u[cbind(bu, 1:ncol(corr_u))])^(-1) *
  #(Tub[cbind(1:B, bu)] - corr_u[cbind(bu, 1:ncol(corr_u))] * Tu)
  corr_u_rows <- corr_u[cbind(bu, bu)]
  taux_2 <- (1 - corr_u_rows)^(-1) * (Tub[cbind(1:B, bu)] - corr_u_rows * Tu)

  tB_u2 <- matrix(1e10, nrow = B, ncol = k_u)
  mask2 <- (1 > corr_u[bu, , drop = FALSE] + tol_r)
  tB_u2[mask2] <- taux_2[mask2]
  th_2_2 <- apply(cbind(tLFu, tB_u2), 1, min) * (Tl < Tu)

  # ----- 5. Combine results -----
  th_1 <- th_1_1 + th_1_2
  th_2 <- th_2_1 + th_2_2

  return(list(th_1 = th_1, th_2 = th_2))
}
