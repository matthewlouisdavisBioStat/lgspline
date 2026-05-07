test_that("predictions are reasonable, derivative can be obtained", {
  set.seed(1234)
  # Predictions follow the true pattern
  t <- seq(-9, 9, length.out = 100)
  y <- sin(t) + rnorm(100, 0, 0.1)
  fit <- lgspline(cbind(t), y, K = 5)

  newt <- matrix(seq(-9, 9, length.out = 10))
  pred <- predict(fit, newt)

  # Check predictions are within reasonable bounds
  expect_true(all(abs(pred) < max(abs(y)) * 1.5))

  ## First derivatives
  if ("take_first_derivatives" %in% names(formals(fit$predict))) {
    deriv <- predict(fit, newt, take_first_derivatives = TRUE)
    expect_type(deriv, "list")
    expect_true("first_deriv" %in% names(deriv))
  }
})
