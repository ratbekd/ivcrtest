test_that("iv_cr_test runs end-to-end", {
  set.seed(123)
  n <- 300

  # Simulate a valid IV structure
  z <- rnorm(n)
  h1 <- rnorm(n)
  h2 <- rnorm(n)
  x <- 0.7 * z + 0.3 * h1 + 0.2 * h2 + rnorm(n)
  y <- 1.5 * x + 0.5 * h1 + 0.3 * h2 + rnorm(n)

  data <- data.frame(y = y, x = x, z = z, h1 = h1, h2 = h2)

  res <- iv_cr_test(
    data = data,
    alpha = 0.05,
    seed = 123,
    rxu_range = c(-0.8,0),
    X = "x", Y = "y", Z = "z",
    H = c("h1", "h2"),
    n = n,
    k = -1
  )
  #
  expect_s3_class(res, "data.frame")
  expect_true(all(c( "CI_Bei", "CI_simple") %in% names(res)))
})
