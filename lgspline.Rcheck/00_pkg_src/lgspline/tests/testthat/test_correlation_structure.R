test_that("basic correlation structure runs without error", {
  set.seed(1234)
  x <- seq(-9, 9, length.out = 1000)
  y <- sin(x) + rnorm(1000, 0, 0.1)
  model_fit <- lgspline(x, y, K = 1, standardize_response = FALSE,
                  VhalfInv = diag(length(y)))
  expect_error(model_fit$model_fit$VhalfInv_params_estimates, NA)
})


test_that("built-in Gaussian correlation structures provide Vhalf_fxn", {
  set.seed(123)

  id <- rep(1:12, each = 4)
  x <- rep(seq(0, 1, length.out = 4), 12)
  y <- sin(2 * pi * x) + rnorm(length(x), 0, 0.05)

  fit <- lgspline(cbind(x), y,
                  K = 1,
                  correlation_id = id,
                  correlation_structure = "exchangeable",
                  include_warnings = FALSE)

  fit_core <- if(!is.null(fit$model_fit)) fit$model_fit else fit

  expect_true(is.function(fit_core$Vhalf_fxn))
  expect_false(is.null(fit_core$VhalfInv_params_estimates))

  Vhalf <- fit_core$Vhalf_fxn(fit_core$VhalfInv_params_estimates)
  expect_true(is.matrix(Vhalf))
  expect_equal(dim(Vhalf), c(length(y), length(y)))
  expect_true(all(is.finite(Vhalf)))
})

test_that("exchangeable correlation supports unconstrained K = 0 GLM fits", {
  set.seed(456)

  id <- rep(1:18, each = 4)
  x <- rep(seq(0, 1, length.out = 4), 18)
  eta <- 0.25 + 0.9 * x
  y <- rpois(length(x), lambda = exp(eta))

  expect_error({
    fit <- lgspline(cbind(x), y,
                    family = quasipoisson(),
                    K = 0,
                    opt = FALSE,
                    correlation_id = id,
                    correlation_structure = "exchangeable",
                    include_warnings = FALSE)
  }, NA)

  fit_core <- if(!is.null(fit$model_fit)) fit$model_fit else fit

  expect_true(is.finite(fit_core$sigmasq_tilde))
  expect_true(all(is.finite(fit_core$varcovmat)))
})

test_that("weighted correlated Gaussian varcov and logLik use the whitened-system weights once", {
  set.seed(20260314)

  n <- 10
  x <- seq(-1, 1, length.out = n)
  y <- 1 + 0.4 * x + rnorm(n, 0, 0.03)
  w <- seq(0.8, 1.7, length.out = n)
  rho <- 0.35
  V <- rho ^ abs(outer(seq_len(n), seq_len(n), `-`))
  VhalfInv <- t(chol(solve(V)))

  fit <- lgspline(
    cbind(x),
    y,
    K = 0,
    opt = FALSE,
    VhalfInv = VhalfInv,
    observation_weights = w,
    just_linear_without_interactions = 1,
    standardize_response = FALSE,
    include_warnings = FALSE
  )

  fit_core <- if(!is.null(fit$model_fit)) fit$model_fit else fit

  X_full <- collapse_block_diagonal(fit_core$X)[unlist(fit_core$og_order), , drop = FALSE]
  VinvhalfX <- fit_core$VhalfInv %*% X_full
  VinvhalfX <- t(t(VinvhalfX) * sqrt(fit_core$weights))

  Lambda_full <- collapse_block_diagonal(
    lapply(seq_len(fit_core$K + 1), function(k) {
      if(length(fit_core$penalties$L_partition_list) == (fit_core$K + 1)) {
        fit_core$penalties$Lambda + fit_core$penalties$L_partition_list[[k]]
      } else {
        fit_core$penalties$Lambda
      }
    })
  )

  gram_expected <- crossprod(VinvhalfX) + Lambda_full
  varcov_expected <- solve(gram_expected) * fit_core$sigmasq_tilde

  resid_w <- c(fit_core$VhalfInv %*% cbind(fit_core$y - fit_core$ytilde))
  logdet_VhalfInv <- determinant(fit_core$VhalfInv, logarithm = TRUE)$modulus[[1L]]
  ll_expected <- -0.5 * fit_core$N * log(2 * pi * fit_core$sigmasq_tilde) +
    logdet_VhalfInv -
    0.5 * sum(fit_core$weights * resid_w^2) / fit_core$sigmasq_tilde

  expect_equal(fit_core$varcovmat, varcov_expected, tolerance = 1e-4)
  expect_equal(as.numeric(logLik(fit_core, include_prior = FALSE)), ll_expected,
               tolerance = 1e-4)
})
