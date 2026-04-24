test_that("lgspline handles basic GLM and quadratic programming constraints", {
  set.seed(1234)

  ## Generate test data
  t <- seq(-9, 9, length.out = 250) # Reduced size for testing

  ## Helper functions
  slinky <- function(t) {
    (50 * cos(t * 2) + -2 * t^2 + (0.25 * t)^4 + 80)
  }

  coil <- function(t) {
    (100 * cos(t * 2) + -1.5 * t^2 + (0.1 * t)^4 + (0.05 * t^3) +
       (-0.01 * t^5) + (0.00002 * t^6) - (0.000001 * t^7) + 100)
  }

  exponential_log <- function(t) {
    unlist(sapply(t, function(tt) {
      if (tt <= 1) {
        100 * (exp(tt) - exp(1))
      } else {
        100 * (log(tt))
      }
    }))
  }

  ## Combined function
  fxn <- function(t) {
    slinky(t) + coil(t) + 0.5 * exponential_log(t) + 25 * t
  }

  ## Create dataset
  dat <- cbind(t, NA)
  colnames(dat) <- c('t', 'y')

  ## Transform mean for quasi-poisson
  m <- mean(fxn(t))
  s <- sd(fxn(t))
  mu <- (fxn(t) - m) / s
  mu <- mu + abs(min(mu)) + 1

  ## Generate poisson responses
  set.seed(1234)  # For reproducibility
  dat[,'y'] <- sapply(mu, function(m) rpois(1, m))

  ## Test fitting with monotonicity constraint and some non-default settings
  fit <- lgspline(dat[,'t', drop=FALSE],
                  dat[,'y'],
                  K = 1,
                  opt = FALSE,
                  wiggle_penalty = 1e-2,
                  flat_ridge_penalty = 1e-2,
                  unique_penalty_per_partition = FALSE,
                  tol = 1e-2,
                  qp_range_lower = 1,
                  qp_monotonic_increase = TRUE,
                  family = quasipoisson())

  ## Basic checks
  expect_s3_class(fit, "lgspline")
  expect_length(fit$B, fit$K + 1)
  expect_true(all(!is.na(fit$ytilde)))

  ## Check monotonicity constraint
  newt <- matrix(sort(t))
  preds <- fit$predict(newt)
  diffs <- diff(preds)
  expect_true(all(diffs >= -1e-2))  # Allow for numerical imprecision

  ## Check range constraint
  expect_true(all(preds > 1-1e-2)) # Allow for numeric imprecision

  ## Test GLM-specific components
  expect_equal(fit$family$family, quasipoisson()$family)
  expect_equal(fit$family$link, quasipoisson()$link)

  ## Test plotting
  expect_error(plot(fit, show_formulas = TRUE,
                    text_size_formula = 2), NA)

  ## Test predictions
  expect_length(predict(fit, matrix(t[1:10])), 10)
})

test_that("Basic lgspline handles logistic regression without constraints", {
  t <- seq(-3, 3, length.out = 250)

  ## Binary response (logistic regression)
  y_bin <- rbinom(250, 1, plogis(sin(t)))
  fit_bin <- lgspline(cbind(t),
                      unique_penalty_per_partition = FALSE,
                      log_initial_flat = 1,
                      log_initial_wiggle = 1e-1,
                      y_bin,
                      K = 10,
                      opt = FALSE,
                      iterate_tune = FALSE,
                      iterate_final_fit = FALSE,
                      include_constrain_first_deriv = FALSE,
                      include_constrain_second_deriv = FALSE,
                      include_constrain_fitted = FALSE,
                      family = quasibinomial())

  ## Class is right
  expect_s3_class(fit_bin, "lgspline")

  ## Family is right
  expect_equal(fit_bin$family$family, quasibinomial()$family)


  ## Functionality works
  generate_posterior(fit_bin, draw_dispersion = FALSE)
  print(summary(fit_bin))
  fit_bin$find_extremum(minimize = TRUE)

  ## Range is right
  preds <- predict(fit_bin, new_predictors = cbind(sample(t)+rnorm(length(t),
                                                                   0,
                                                                   0.00001)))
  expect_true(all(abs(preds) <= 1))
})

test_that("lgspline handles various quadratic programming constraints", {
  t <- seq(-3, 3, length.out = 100)
  y <- exp(t) + rnorm(100, 0, 0.1)

  ## Test monotone increasing
  fit_inc <- lgspline(cbind(t), y,
                      K = 2,
                      opt = FALSE,
                      qp_monotonic_increase = TRUE)
  preds_inc <- predict(fit_inc, matrix(sort(t)))
  expect_true(all(diff(preds_inc) >= -1e-10))

  ## Test monotone decreasing
  y_dec <- -y
  fit_dec <- lgspline(cbind(t), y_dec,
                      K = 0,
                      opt = FALSE,
                      qp_monotonic_decrease = TRUE)
  preds_dec <- predict(fit_dec, matrix(sort(t)))
  expect_true(all(diff(preds_dec) <= 1e-10))

  ## Test bounded range
  y_bound <- y - mean(y)
  fit_bound <- lgspline(cbind(t), y_bound,
                        K = 1,
                        qp_range_lower = -1,
                        qp_range_upper = 1,
                        opt = FALSE)
  preds_bound <- predict(fit_bound, matrix(t))
  expect_true(all(preds_bound >= -1.01)) # Allow small numerical error
  expect_true(all(preds_bound <= 1.01))
})

test_that("Gaussian QP constraints are invariant to response standardization", {
  set.seed(20260423)

  t <- seq(-4, 4, length.out = 90)
  y <- 2 + exp(-0.4 * t) + 0.15 * t^2 + rnorm(length(t), 0, 0.15)

  args <- list(
    predictors = cbind(t),
    y = y,
    K = 2,
    opt = FALSE,
    qp_range_lower = 0,
    qp_negative_derivative = TRUE,
    qp_positive_2ndderivative = TRUE,
    include_cubic_terms = FALSE,
    include_quartic_terms = FALSE,
    include_constrain_second_deriv = FALSE,
    include_warnings = FALSE
  )

  fit_std <- do.call(lgspline, c(args, list(standardize_response = TRUE)))
  fit_raw <- do.call(lgspline, c(args, list(standardize_response = FALSE)))

  pred_std <- predict(fit_std, cbind(t))
  pred_raw <- predict(fit_raw, cbind(t))

  expect_equal(pred_std, pred_raw, tolerance = 1e-6)
  expect_true(all(pred_std >= -1e-8))

  deriv_std <- predict(fit_std, cbind(t), take_first_derivatives = TRUE)
  deriv_raw <- predict(fit_raw, cbind(t), take_first_derivatives = TRUE)

  expect_equal(unlist(deriv_std$first_deriv[[1]]),
               unlist(deriv_raw$first_deriv[[1]]),
               tolerance = 1e-6)
})
