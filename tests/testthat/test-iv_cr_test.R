# test-iv_cr_test.R

library(testthat)

test_that("iv_cr_test returns expected CI structure", {
  set.seed(123)

  # Simulate valid data for IV regression test
  n <- 300
  z <- rnorm(n)
  h1 <- rnorm(n)
  h2 <- rnorm(n)
  h3 <- rnorm(n)

  # Endogenous regressor
  x <- 0.7 * z + 0.4 * h1 + 0.3 * h2 + 0.2 * h3 + rnorm(n, sd = 0.5)

  # Outcome depending on x and exogenous H
  y <- 1.5 * x + 0.3 * h1 + 0.2 * h2 + 0.1 * h3 + rnorm(n)

  df <- data.frame(y = y, x = x, z = z, h1 = h1, h2 = h2, h3 = h3)

  # Variable names (strings)
  Y <- "y"
  X <- "x"
  Z <- "z"
  H <- c("h1", "h2", "h3")

  # Run the function
  res <- iv_cr_test(data = df, X = X, Y = Y, H = H, Z = Z, n = nrow(df), k = -1, alpha = 0.05,
                    seed = 123,
                    rxu_range = c(0, 0.8),
                    bias_mc = FALSE,
                    mc_B = 500)
    # --- Expectations ---
  # Check that result is a data frame
  expect_s3_class(res, "data.frame")

  # Check that it has required columns
  expected_cols <- c("Z", "plug_in", "CI_Bei", "CI_simple",
                     "Zero_in_CI","p_zero","n")

  expect_true(all(expected_cols %in% colnames(res)))

  # Check result has 1 row (one IV test)
  expect_equal(nrow(res), 1)

  # Check that confidence intervals are strings formatted like [a,b]
  expect_true(grepl("^\\[.*\\]$", res$CI_Bei))
  expect_true(grepl("^\\[.*\\]$", res$CI_simple))
  # expect_true(grepl("^\\[.*\\]$", res$p_zero))
  expect_true(is.numeric(res$p_zero) &&
                !is.na(res$p_zero) &&
                res$p_zero >= 0 &&
                res$p_zero <= 1)
})
