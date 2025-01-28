test_that("predictions are reasonable, derivative can be obtained", {
  set.seed(1234)
  # Test predictions follow true function pattern
  x <- seq(-9, 9, length.out = 100)
  y <- sin(x) + rnorm(100, 0, 0.1)
  fit <- lgspline(cbind(x), y, K = 5)

  newx <- matrix(seq(-9, 9, length.out = 10))
  pred <- predict(fit, newx)

  # Check predictions are within reasonable bounds
  expect_true(all(abs(pred) < max(abs(y)) * 1.5))

  ## Test first derivatives
  if ("take_first_derivatives" %in% names(formals(fit$predict))) {
    deriv <- predict(fit, newx, take_first_derivatives = TRUE)
    expect_type(deriv, "list")
    expect_true("first_deriv" %in% names(deriv))
  }
})

