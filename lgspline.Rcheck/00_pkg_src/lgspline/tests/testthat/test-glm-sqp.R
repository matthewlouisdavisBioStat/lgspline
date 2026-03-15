test_that("lgspline handles basic GLM and quadratic programming constraints", {
  set.seed(1234)

  ## Generate test data
  x <- seq(-9, 9, length.out = 250) # Reduced size for testing

  ## Helper functions
  slinky <- function(x) {
    (50 * cos(x * 2) + -2 * x^2 + (0.25 * x)^4 + 80)
  }

  coil <- function(x) {
    (100 * cos(x * 2) + -1.5 * x^2 + (0.1 * x)^4 + (0.05 * x^3) +
       (-0.01 * x^5) + (0.00002 * x^6) - (0.000001 * x^7) + 100)
  }

  exponential_log <- function(x) {
    unlist(sapply(x, function(xx) {
      if (xx <= 1) {
        100 * (exp(xx) - exp(1))
      } else {
        100 * (log(xx))
      }
    }))
  }

  ## Combined function
  fxn <- function(x) {
    slinky(x) + coil(x) + 0.5 * exponential_log(x) + 25 * x
  }

  ## Create dataset
  dat <- cbind(x, NA)
  colnames(dat) <- c('x', 'y')

  ## Transform mean for quasi-poisson
  m <- mean(fxn(x))
  s <- sd(fxn(x))
  mu <- (fxn(x) - m) / s
  mu <- mu + abs(min(mu)) + 1

  ## Generate poisson responses
  set.seed(1234)  # For reproducibility
  dat[,'y'] <- sapply(mu, function(m) rpois(1, m))

  ## Test fitting with monotonicity constraint and some non-default settings
  fit <- lgspline(dat[,'x', drop=FALSE],
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
  newx <- matrix(sort(x))
  preds <- fit$predict(newx)
  diffs <- diff(preds)
  expect_true(all(diffs >= -1e-10))  # Allow for numerical imprecision

  ## Check range constraint
  expect_true(all(preds > 1-1e-10)) # Allow for numeric imprecision

  ## Test GLM-specific components
  expect_equal(fit$family$family, quasipoisson()$family)
  expect_equal(fit$family$link, quasipoisson()$link)

  ## Test plotting
  expect_error(plot(fit, show_formulas = TRUE,
                    text_size_formula = 2), NA)

  ## Test predictions
  expect_length(predict(fit, matrix(x[1:10])), 10)
})

test_that("Basic lgspline handles logistic regression without constraints", {
  x <- seq(-3, 3, length.out = 250)

  ## Binary response (logistic regression)
  y_bin <- rbinom(250, 1, plogis(sin(x)))
  fit_bin <- lgspline(cbind(x),
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
  preds <- predict(fit_bin, new_predictors = cbind(sample(x)+rnorm(length(x),
                                                             0,
                                                             0.00001)))
  expect_true(all(abs(preds) <= 1))
})

test_that("lgspline handles various quadratic programming constraints", {
  x <- seq(-3, 3, length.out = 100)
  y <- exp(x) + rnorm(100, 0, 0.1)

  ## Test monotone increasing
  fit_inc <- lgspline(cbind(x), y,
                      K = 2,
                      opt = FALSE,
                      qp_monotonic_increase = TRUE)
  preds_inc <- predict(fit_inc, matrix(sort(x)))
  expect_true(all(diff(preds_inc) >= -1e-10))

  ## Test monotone decreasing
  y_dec <- -y
  fit_dec <- lgspline(cbind(x), y_dec,
                      K = 0,
                      opt = FALSE,
                      qp_monotonic_decrease = TRUE)
  preds_dec <- predict(fit_dec, matrix(sort(x)))
  expect_true(all(diff(preds_dec) <= 1e-10))

  ## Test bounded range
  y_bound <- y - mean(y)
  fit_bound <- lgspline(cbind(x), y_bound,
                       K = 1,
                       qp_range_lower = -1,
                       qp_range_upper = 1,
                       opt = FALSE)
  preds_bound <- predict(fit_bound, matrix(x))
  expect_true(all(preds_bound >= -1.01)) # Allow small numerical error
  expect_true(all(preds_bound <= 1.01))
})
