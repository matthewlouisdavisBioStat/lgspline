test_that("process_input resolves names, offsets, and encoded factor groups", {
  df <- data.frame(
    y = c(1, 2, 3, 4, 5, 6),
    t = seq(-1, 1, length.out = 6),
    z = c(10, 11, 12, 13, 14, 15),
    grp = factor(c("a", "b", "c", "a", "b", "c"))
  )

  processed <- process_input(
    predictors = y ~ spl(t) + grp + offset(z),
    data = df,
    just_linear_without_interactions = "z",
    do_not_cluster_on_these = "z",
    auto_encode_factors = TRUE,
    include_warnings = FALSE
  )

  expect_true(is.matrix(processed$predictors))
  expect_equal(length(processed$offset), 1)
  expect_equal(length(processed$do_not_cluster_on_these), 1)
  expect_identical(processed$offset, processed$do_not_cluster_on_these)
  expect_true("grp" %in% names(processed$factor_groups))
  expect_true(length(processed$factor_groups$grp) > 1)
  expect_true(all(processed$factor_groups$grp %in%
                    processed$just_linear_without_interactions))
})

test_that("integration and prior helper methods return coherent values", {
  t <- seq(-1, 2, length.out = 100)
  y <- 2 * t + 1

  fit <- lgspline(
    t,
    y,
    K = 1,
    opt = FALSE,
    just_linear_without_interactions = 1,
    standardize_response = FALSE
  )
  fit_linear_only <- lgspline(
    t,
    y,
    K = 1,
    opt = FALSE,
    just_linear_with_interactions = 1,
    standardize_response = FALSE
  )

  integ_resp <- integrate(fit, lower = -1, upper = 2, n_quad = 20)
  integ_link <- integrate(fit, lower = -1, upper = 2, n_quad = 20,
                          link_scale = TRUE)
  integ_linear_only <- integrate(fit_linear_only, lower = -1, upper = 2,
                                 n_quad = 20)

  expect_equal(integ_resp, 6, tolerance = 1e-4)
  expect_equal(integ_link, integ_resp, tolerance = 1e-4)
  expect_equal(integ_linear_only, integ_resp, tolerance = 1e-4)

  lp_with_const <- prior_loglik(fit, include_constant = TRUE)
  lp_no_const <- prior_loglik(fit, include_constant = FALSE)

  part_penalties <- fit$penalties$L_partition_list
  if(length(part_penalties) == 0){
    part_penalties <- lapply(seq_len(fit$K + 1), function(k) 0)
  }
  expected_const <- sum(vapply(seq_len(fit$K + 1), function(k) {
    Lambda_total <- fit$penalties$Lambda + part_penalties[[k]]
    -0.5 * (
      length(fit$B_raw[[k]]) * log(2 * pi) +
        length(fit$B_raw[[k]]) * log(fit$sigmasq_tilde) -
        as.numeric(determinant(Lambda_total, logarithm = TRUE)$modulus)
    )
  }, numeric(1)))

  expect_true(is.numeric(lp_with_const) && length(lp_with_const) == 1)
  expect_true(is.numeric(lp_no_const) && length(lp_no_const) == 1)
  expect_true(is.finite(lp_with_const))
  expect_true(is.finite(lp_no_const))
  expect_equal(lp_with_const - lp_no_const, expected_const, tolerance = 1e-4)
})

test_that("leave_one_out matches explicit Gaussian refits in a linear case", {
  t <- seq(-2, 2, length.out = 8)
  y <- 1 + 0.5 * t + c(0.1, -0.05, 0.02, -0.03, 0.04, -0.02, 0.01, -0.04)

  fit <- lgspline(
    cbind(t),
    y,
    K = 0,
    opt = FALSE,
    just_linear_without_interactions = 1,
    standardize_response = FALSE
  )

  loo_fast <- leave_one_out(fit)

  loo_refit <- sapply(seq_along(y), function(i) {
    refit <- lgspline(
      cbind(t[-i]),
      y[-i],
      K = 0,
      opt = FALSE,
      just_linear_without_interactions = 1,
      standardize_response = FALSE
    )
    as.numeric(predict(refit, cbind(t[i])))
  })

  expect_equal(loo_fast, loo_refit, tolerance = 1e-4)
})
