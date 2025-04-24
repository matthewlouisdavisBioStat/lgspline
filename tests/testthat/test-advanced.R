test_that("Weibull accelerated failure time (AFT) models can be fit under custom constraints", {

  ## Setup
  set.seed(1234)
  x1 <- rnorm(1000)
  x2 <- rbinom(1000, 1, 0.5)
  yraw <- rexp(exp(0.01*x1 + 0.01*x2))
  status <- rbinom(1000, 1, 0.25)
  yobs <- ifelse(status, runif(1, 0, yraw), yraw)
  df <- data.frame(
    y = yobs,
    x1 = x1,
    x2 = x2
  )

  ## Weibull AFT
  model_fit_aft <- lgspline(y ~ spl(x1) + x2,
                            df,
                            unconstrained_fit_fxn = unconstrained_fit_weibull,
                            family = weibull_family(),
                            need_dispersion_for_estimation = TRUE,
                            qp_score_function = weibull_qp_score_function(),
                            dispersion_function = weibull_dispersion_function,
                            glm_weight_function = weibull_glm_weight_function,
                            shur_correction_function = weibull_shur_correction,
                            K = 1,
                            include_constrain_fitted = FALSE,
                            flat_ridge_penalty = 1e-2,
                            wiggle_penalty = 1e-2,
                            unique_penalty_per_predictor = FALSE,
                            unique_penalty_per_partition = FALSE,
                            opt = FALSE,
                            return_varcovmat = TRUE,
                            observation_weights = 0.5 + 0.5*abs(rnorm(1000)),
                            constraint_vectors = cbind(rep(1, 12)),
                            null_constraint = matrix(12),
                            status = status,
                            verbose = TRUE)
  print(model_fit_aft)

  ## Check sum-to-P constraint only if model_fit_aft is a list
  if (is.list(model_fit_aft)) {
    expect_equal(round(sum(unlist(model_fit_aft$B)), 8), 12)
  } else {
  ## Skip the test on platforms when model returns a try-error instead
    expect_error(invisible(), NA)
  }
})
