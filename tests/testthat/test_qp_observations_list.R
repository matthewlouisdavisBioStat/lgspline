fit_core <- function(fit) {
  if (!is.null(fit$model_fit)) fit$model_fit else fit
}


build_qp <- function(fit_core, qp_observations = NULL, ...) {
  lgspline:::process_qp(
    X = fit_core$X,
    K = fit_core$K,
    p_expansions = fit_core$p,
    order_list = fit_core$order_list,
    colnm_expansions = fit_core$raw_expansion_names,
    expansion_scales = fit_core$expansion_scales,
    power1_cols = fit_core$power1_cols,
    power2_cols = fit_core$power2_cols,
    nonspline_cols = fit_core$nonspline_cols,
    interaction_single_cols = fit_core$interaction_single_cols,
    interaction_quad_cols = fit_core$interaction_quad_cols,
    triplet_cols = fit_core$triplet_cols,
    include_2way_interactions =
      fit_core$.fit_call_args$include_2way_interactions,
    include_3way_interactions =
      fit_core$.fit_call_args$include_3way_interactions,
    include_quadratic_interactions =
      fit_core$.fit_call_args$include_quadratic_interactions,
    family = fit_core$family,
    mean_y = fit_core$mean_y,
    sd_y = fit_core$sd_y,
    N_obs = fit_core$N,
    qp_observations = qp_observations,
    og_cols = colnames(fit_core$.fit_call_args$predictors),
    include_warnings = FALSE,
    ...
  )
}


expect_qp_equal <- function(qp1, qp2) {
  expect_equal(qp1$qp_Amat, qp2$qp_Amat)
  expect_equal(qp1$qp_bvec, qp2$qp_bvec)
  expect_equal(qp1$qp_meq, qp2$qp_meq)
}


test_that("keyed qp_observations reproduces the legacy shared subset path", {
  set.seed(1234)

  Time <- seq(-4, 4, length.out = 80)
  y <- exp(-0.35 * Time) + rnorm(length(Time), 0, 0.02)
  x <- cbind(Time = Time)
  obs <- which(Time > 0)

  args <- list(
    predictors = x,
    y = y,
    K = 2,
    opt = FALSE,
    standardize_response = FALSE,
    qp_range_lower = 0,
    qp_negative_derivative = "Time",
    qp_positive_2ndderivative = "Time",
    include_cubic_terms = FALSE,
    include_quartic_terms = FALSE,
    include_constrain_second_deriv = FALSE,
    include_warnings = FALSE
  )

  set.seed(1234)
  fit_vec <- do.call(lgspline::lgspline, c(args, list(
    qp_observations = obs
  )))
  set.seed(1234)
  fit_key <- do.call(lgspline::lgspline, c(args, list(
    qp_observations = list(
      "qp_range_lower" = obs,
      "Time:qp_negative_derivative" = obs,
      "Time:qp_positive_2ndderivative" = obs
    )
  )))

  fit_vec_core <- fit_core(fit_vec)
  fit_key_core <- fit_core(fit_key)
  qp_vec <- build_qp(
    fit_vec_core,
    qp_observations = obs,
    qp_range_lower = 0,
    qp_negative_derivative = "Time",
    qp_positive_2ndderivative = "Time"
  )
  qp_key <- build_qp(
    fit_key_core,
    qp_observations = list(
      "qp_range_lower" = obs,
      "Time:qp_negative_derivative" = obs,
      "Time:qp_positive_2ndderivative" = obs
    ),
    qp_range_lower = 0,
    qp_negative_derivative = "Time",
    qp_positive_2ndderivative = "Time"
  )

  expect_qp_equal(qp_vec, qp_key)
  expect_equal(stats::predict(fit_vec, x),
               stats::predict(fit_key, x),
               tolerance = 1e-8)
})


test_that("keyed qp_observations combines per-constraint subsets correctly", {
  set.seed(1234)

  Time <- seq(-2, 6, length.out = 90)
  y <- exp(0.25 * Time) + rnorm(length(Time), 0, 0.02)
  x <- cbind(Time = Time)

  obs_range <- which(Time >= 0)
  obs_pos1 <- which(Time > 1)
  obs_pos2 <- which(Time > 2)

  fit <- lgspline::lgspline(
    x,
    y,
    K = 2,
    opt = FALSE,
    standardize_response = FALSE,
    qp_observations = list(
      "qp_range_lower" = obs_range,
      "Time:qp_positive_derivative" = obs_pos1,
      "Time:qp_positive_2ndderivative" = obs_pos2
    ),
    qp_range_lower = 0,
    qp_positive_derivative = "Time",
    qp_positive_2ndderivative = "Time",
    include_cubic_terms = FALSE,
    include_quartic_terms = FALSE,
    include_constrain_second_deriv = FALSE,
    include_warnings = FALSE
  )

  fit_core_obj <- fit_core(fit)
  qp_key <- build_qp(
    fit_core_obj,
    qp_observations = list(
      "qp_range_lower" = obs_range,
      "Time:qp_positive_derivative" = obs_pos1,
      "Time:qp_positive_2ndderivative" = obs_pos2
    ),
    qp_range_lower = 0,
    qp_positive_derivative = "Time",
    qp_positive_2ndderivative = "Time"
  )
  qp_range <- build_qp(
    fit_core_obj,
    qp_observations = obs_range,
    qp_range_lower = 0
  )
  qp_pos1 <- build_qp(
    fit_core_obj,
    qp_observations = obs_pos1,
    qp_positive_derivative = "Time"
  )
  qp_pos2 <- build_qp(
    fit_core_obj,
    qp_observations = obs_pos2,
    qp_positive_2ndderivative = "Time"
  )

  expect_equal(
    qp_key$qp_Amat,
    cbind(qp_range$qp_Amat, qp_pos1$qp_Amat, qp_pos2$qp_Amat)
  )
  expect_equal(
    qp_key$qp_bvec,
    c(qp_range$qp_bvec, qp_pos1$qp_bvec, qp_pos2$qp_bvec)
  )
  expect_equal(qp_key$qp_meq, 0)
})


test_that("keyed qp_observations can vary by variable and ignore unknown keys", {
  set.seed(1234)

  n <- 110
  Time <- sort(runif(n, -2, 4))
  Dose <- sort(runif(n, 0.5, 2.5))
  y <- exp(0.35 * Time) + Dose + rnorm(n, 0, 0.02)
  x <- cbind(Time = Time, Dose = Dose)

  obs_time <- which(Time > 0)
  obs_dose <- which(Dose < 1.8)

  expect_warning(
    fit <- lgspline::lgspline(
      x,
      y,
      K = 1,
      opt = FALSE,
      standardize_response = FALSE,
      qp_observations = list(
        "Time:qp_positive_derivative" = obs_time,
        "Dose:qp_positive_derivative" = obs_dose,
        "Time:qp_banana" = obs_time
      ),
      qp_positive_derivative = c("Time", "Dose"),
      include_cubic_terms = FALSE,
      include_quartic_terms = FALSE,
      include_constrain_second_deriv = FALSE
    ),
    "not a known constraint type"
  )

  fit_core_obj <- fit_core(fit)
  qp_key <- build_qp(
    fit_core_obj,
    qp_observations = list(
      "Time:qp_positive_derivative" = obs_time,
      "Dose:qp_positive_derivative" = obs_dose
    ),
    qp_positive_derivative = c("Time", "Dose")
  )
  qp_time <- build_qp(
    fit_core_obj,
    qp_observations = obs_time,
    qp_positive_derivative = "Time"
  )
  qp_dose <- build_qp(
    fit_core_obj,
    qp_observations = obs_dose,
    qp_positive_derivative = "Dose"
  )

  expect_equal(
    qp_key$qp_Amat,
    cbind(qp_time$qp_Amat, qp_dose$qp_Amat)
  )
  expect_equal(
    qp_key$qp_bvec,
    c(qp_time$qp_bvec, qp_dose$qp_bvec)
  )
})
