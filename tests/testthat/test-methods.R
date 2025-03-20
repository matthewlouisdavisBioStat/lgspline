test_that("basic S3 methods work correctly", {
  x <- seq(-9, 9, length.out = 100)
  y <- sin(x) + rnorm(100, 0, 0.1) + 0.1*x^2
  fit <- lgspline(cbind(x), y, K = 5)

  ## Test print method
  expect_output(print(fit), "Lagrangian Multiplier Smoothing Spline Model")

  ## Test predict method
  newx <- matrix(seq(-9, 9, length.out = 10))
  pred <- predict(fit, newx)
  expect_length(pred, 10)

  ## Test coef method
  coefs <- coef(fit)
  expect_type(coefs, "list")
  expect_length(coefs, fit$K + 1)

  ## Test find extremum
  extr <- find_extremum(fit)
  expect_type(extr, "list")
  expect_length(extr, 2)

  ## Test generate posterior
  extr <- generate_posterior(fit)
  expect_type(extr, "list")
  expect_length(extr, 2)

})
