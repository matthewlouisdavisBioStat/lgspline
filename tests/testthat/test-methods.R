test_that("basic S3 methods work correctly", {
  t <- seq(-9, 9, length.out = 100)
  y <- sin(t) + rnorm(100, 0, 0.1) + 0.1*t^2
  fit <- lgspline(cbind(t), y, K = 5)
  ## Test print method
  expect_output(print(fit), "Lagrangian Multiplier Smoothing Spline Model")
  ## Test predict method
  newt <- matrix(seq(-9, 9, length.out = 10))
  pred <- predict(fit, newt)
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
test_that("generate_posterior forwards enforce_qp_constraints to the stored sampler", {
  set.seed(2026)
  t <- seq(0, 1, length.out = 40)
  y <- exp(t) + rnorm(length(t), 0, 0.02)
  fit <- lgspline(cbind(t), y,
                  K = 1,
                  opt = FALSE,
                  qp_monotonic_increase = TRUE,
                  include_warnings = FALSE)
  fit_core <- if(!is.null(fit$model_fit)) fit$model_fit else fit
  sampler <- fit_core$generate_posterior
  expect_true("enforce_qp_constraints" %in% names(formals(sampler)))
  set.seed(99)
  draw_direct <- sampler(draw_dispersion = FALSE,
                         num_draws = 1,
                         enforce_constraints = TRUE,
                         max_slice_iterations = 200)
  set.seed(99)
  draw_public <- generate_posterior(fit,
                                    draw_dispersion = FALSE,
                                    num_draws = 1,
                                    enforce_qp_constraints = TRUE,
                                    max_slice_iterations = 200)
  expect_equal(draw_public$post_draw_coefficients,
               draw_direct$post_draw_coefficients)
})
