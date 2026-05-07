local_null_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

build_tuning_env_from_fit <- function(fit,
                                      tuning_criterion = "loo",
                                      gcv_gamma = 1.4) {
  X_std <- lapply(fit$X, fit$std_X)
  smoothing_spline_penalty <- lgspline:::get_2ndDerivPenalty(
    colnm_expansions = fit$raw_expansion_names,
    C = fit$X[[1]],
    power1_cols = fit$power1_cols,
    power2_cols = fit$power2_cols,
    power3_cols = fit$power3_cols,
    power4_cols = fit$power4_cols,
    interaction_single_cols = fit$interaction_single_cols,
    interaction_quad_cols = fit$interaction_quad_cols,
    triplet_cols = fit$triplet_cols,
    p_expansions = fit$p
  )
  y_vec <- fit$y
  y_list <- lapply(fit$order_list, function(idx) y_vec[idx])
  X_gram <- lapply(X_std, crossprod)
  Xy <- lapply(seq_along(X_std), function(k) crossprod(X_std[[k]], cbind(y_list[[k]])))
  delta <- fit$.fit_call_args$delta
  if (is.null(delta) && paste0(fit$family)[2] != "identity") {
    delta <- lgspline:::.compute_tuning_delta(
      family = fit$family,
      unl_y = y_vec,
      N_obs = fit$N,
      observation_weights = NULL,
      opt = TRUE
    )
  }
  if (is.null(delta)) {
    delta <- 0
  }

  A <- fit$A
  if (is.null(A)) {
    A <- cbind(rep(0, (fit$K + 1) * fit$p))
    A <- cbind(A, A)
  }
  R_constraints <- ncol(A)

  observation_weights <- fit$weights
  if (is.null(observation_weights)) {
    observation_weights <- rep(1, fit$N)
  }
  observation_weights_list <- lapply(fit$order_list, function(idx) {
    observation_weights[idx]
  })
  homogenous_weights <- length(unique(observation_weights)) == 1

  flat_cols <- c()
  flat_terms <- local_null_coalesce(
    fit$.fit_call_args$just_linear_without_interactions,
    c()
  )
  if (length(flat_terms) > 0) {
    flat_cols <- which(fit$raw_expansion_names %in% paste0("_", flat_terms, "_"))
  }
  use_blockfit <- isTRUE(fit$.fit_call_args$blockfit) &&
    length(flat_cols) > 0 && fit$K > 0

  lgspline:::.build_tuning_env(
    y = y_list,
    X = X_std,
    X_gram = X_gram,
    Xy = Xy,
    smoothing_spline_penalty = smoothing_spline_penalty,
    A = A,
    R_constraints = R_constraints,
    K = fit$K,
    p_expansions = fit$p,
    N_obs = fit$N,
    custom_penalty_mat = fit$.fit_call_args$custom_penalty_mat,
    colnm_expansions = fit$raw_expansion_names,
    unique_penalty_per_predictor = fit$.fit_call_args$unique_penalty_per_predictor,
    unique_penalty_per_partition = fit$.fit_call_args$unique_penalty_per_partition,
    meta_penalty = fit$.fit_call_args$meta_penalty,
    family = fit$family,
    delta = delta,
    order_list = fit$order_list,
    observation_weights = observation_weights_list,
    homogenous_weights = homogenous_weights,
    parallel = FALSE,
    parallel_eigen = FALSE,
    parallel_trace = FALSE,
    parallel_aga = FALSE,
    parallel_matmult = FALSE,
    parallel_unconstrained = FALSE,
    cl = NULL,
    chunk_size = 1,
    num_chunks = 0,
    rem_chunks = 0,
    unconstrained_fit_fxn = fit$.fit_call_args$unconstrained_fit_fxn,
    keep_weighted_Lambda = fit$.fit_call_args$keep_weighted_Lambda,
    iterate = fit$.fit_call_args$iterate_tune,
    qp_score_function = fit$.fit_call_args$qp_score_function,
    quadprog = FALSE,
    qp_Amat = NULL,
    qp_bvec = NULL,
    qp_meq = 0,
    tol = fit$.fit_call_args$tol,
    sd_y = fit$sd_y,
    tuning_criterion = tuning_criterion,
    gcv_gamma = gcv_gamma,
    constraint_value_vectors = local_null_coalesce(fit$constraint_values, cbind()),
    glm_weight_function = fit$.fit_call_args$glm_weight_function,
    schur_correction_function = fit$.fit_call_args$schur_correction_function,
    need_dispersion_for_estimation = fit$.fit_call_args$need_dispersion_for_estimation,
    dispersion_function = fit$.fit_call_args$dispersion_function,
    blockfit = fit$.fit_call_args$blockfit,
    just_linear_without_interactions = flat_terms,
    Vhalf = fit$Vhalf,
    VhalfInv = fit$VhalfInv,
    verbose = FALSE,
    include_warnings = FALSE,
    flat_cols = flat_cols,
    use_blockfit = use_blockfit
  )
}

make_gaussian_linear_fit <- function(wiggle_penalty = 1e-3,
                                     flat_ridge_penalty = 0.2) {
  t <- seq(-2, 2, length.out = 8)
  y <- 1 + 0.5 * t + c(0.1, -0.05, 0.02, -0.03, 0.04, -0.02, 0.01, -0.04)

  fit <- lgspline(
    cbind(t),
    y,
    K = 0,
    opt = FALSE,
    wiggle_penalty = wiggle_penalty,
    flat_ridge_penalty = flat_ridge_penalty,
    unique_penalty_per_predictor = FALSE,
    unique_penalty_per_partition = FALSE,
    just_linear_without_interactions = 1,
    standardize_response = FALSE,
    include_warnings = FALSE
  )

  list(
    t = t,
    y = y,
    fit = fit,
    wiggle_penalty = wiggle_penalty,
    flat_ridge_penalty = flat_ridge_penalty
  )
}

test_that(".compute_loocv matches explicit Gaussian refits at fixed penalties", {
  case <- make_gaussian_linear_fit()
  env <- build_tuning_env_from_fit(case$fit, tuning_criterion = "loo")
  par <- log(c(case$wiggle_penalty, case$flat_ridge_penalty))

  loo_obj <- lgspline:::.compute_loocv(par, c(), env)

  loo_refit <- sapply(seq_along(case$y), function(i) {
    refit <- lgspline(
      cbind(case$t[-i]),
      case$y[-i],
      K = 0,
      opt = FALSE,
      wiggle_penalty = exp(par[1]),
      flat_ridge_penalty = exp(par[2]),
      unique_penalty_per_predictor = FALSE,
      unique_penalty_per_partition = FALSE,
      just_linear_without_interactions = 1,
      standardize_response = FALSE,
      include_warnings = FALSE
    )
    as.numeric(predict(refit, cbind(case$t[i])))
  })

  loo_refit_value <- mean((case$y - loo_refit)^2)

  expect_equal(loo_obj$LOO_u, loo_refit_value, tolerance = 1e-4)
  expect_equal(loo_obj$criterion_value, loo_obj$LOO_u, tolerance = 1e-10)
})

test_that(".compute_loocv_gradient matches finite differences in a Gaussian case", {
  case <- make_gaussian_linear_fit(wiggle_penalty = 5e-3, flat_ridge_penalty = 0.3)
  env <- build_tuning_env_from_fit(case$fit, tuning_criterion = "loo")
  par <- log(c(5e-3, 0.3))

  analytic <- lgspline:::.compute_loocv_gradient(par, c(), env = env)$gradient[1:2]

  h <- 1e-5
  numeric <- vapply(seq_along(par), function(j) {
    step <- rep(0, length(par))
    step[j] <- h
    up <- lgspline:::.compute_loocv(par + step, c(), env)$criterion_value
    down <- lgspline:::.compute_loocv(par - step, c(), env)$criterion_value
    (up - down) / (2 * h)
  }, numeric(1))

  expect_equal(analytic, numeric, tolerance = 1e-4)
})

test_that("default tuning is LOO and gcv_gamma affects only the GCV objective", {
  t <- seq(-2, 2, length.out = 30)
  y <- sin(t) + 0.1 * t^2
  fit_gcv_case <- lgspline(
    cbind(t), y,
    K = 1,
    opt = FALSE,
    wiggle_penalty = 5e-3,
    flat_ridge_penalty = 0.3,
    unique_penalty_per_predictor = FALSE,
    unique_penalty_per_partition = FALSE,
    standardize_response = FALSE,
    include_warnings = FALSE
  )
  env_loo_1 <- build_tuning_env_from_fit(fit_gcv_case, tuning_criterion = "loo", gcv_gamma = 1)
  env_loo_9 <- build_tuning_env_from_fit(fit_gcv_case, tuning_criterion = "loo", gcv_gamma = 9)
  env_gcv_1 <- build_tuning_env_from_fit(fit_gcv_case, tuning_criterion = "gcv", gcv_gamma = 1)
  env_gcv_9 <- build_tuning_env_from_fit(fit_gcv_case, tuning_criterion = "gcv", gcv_gamma = 9)
  par <- log(c(5e-3, 0.3))

  loo_1 <- lgspline:::.compute_loocv(par, c(), env_loo_1)$criterion_value
  loo_9 <- lgspline:::.compute_loocv(par, c(), env_loo_9)$criterion_value
  gcv_1 <- lgspline:::.compute_gcvu(par, c(), env_gcv_1)$criterion_value
  gcv_9 <- lgspline:::.compute_gcvu(par, c(), env_gcv_9)$criterion_value

  expect_equal(loo_1, loo_9, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(gcv_1, gcv_9, tolerance = 1e-10)))

  t <- seq(-1, 1, length.out = 20)
  y <- sin(t) + 0.1 * cos(3 * t)
  fit_default <- lgspline(
    cbind(t), y,
    K = 1,
    opt = TRUE,
    use_custom_bfgs = FALSE,
    gcv_gamma = 9,
    initial_wiggle = c(1e-4, 5e-4),
    initial_flat = c(0.2, 0.5),
    unique_penalty_per_predictor = FALSE,
    unique_penalty_per_partition = FALSE,
    just_linear_without_interactions = 1,
    standardize_response = FALSE,
    include_warnings = FALSE
  )
  fit_explicit_loo <- lgspline(
    cbind(t), y,
    K = 1,
    opt = TRUE,
    use_custom_bfgs = FALSE,
    tuning_criterion = "loo",
    gcv_gamma = 1,
    initial_wiggle = c(1e-4, 5e-4),
    initial_flat = c(0.2, 0.5),
    unique_penalty_per_predictor = FALSE,
    unique_penalty_per_partition = FALSE,
    just_linear_without_interactions = 1,
    standardize_response = FALSE,
    include_warnings = FALSE
  )

  expect_equal(fit_default$penalties$Lambda,
               fit_explicit_loo$penalties$Lambda,
               tolerance = 1e-8)
})

test_that("LOO and GCV tuning helpers stay finite for a non-Gaussian case", {
  set.seed(2)
  t <- seq(-1, 1, length.out = 40)
  eta <- -0.3 + 1.2 * t
  p <- plogis(eta)
  y <- rbinom(length(t), 1, p)

  fit <- lgspline(
    cbind(t),
    y,
    family = binomial(),
    K = 1,
    opt = FALSE,
    wiggle_penalty = 1e-2,
    flat_ridge_penalty = 0.5,
    unique_penalty_per_predictor = FALSE,
    unique_penalty_per_partition = FALSE,
    standardize_response = FALSE,
    include_warnings = FALSE
  )

  env_loo <- build_tuning_env_from_fit(fit, tuning_criterion = "loo")
  env_gcv <- build_tuning_env_from_fit(fit, tuning_criterion = "gcv", gcv_gamma = 1.4)
  par <- log(c(1e-2, 0.5))

  loo_obj <- lgspline:::.compute_loocv(par, c(), env_loo)
  loo_grad <- lgspline:::.compute_loocv_gradient(par, c(), env = env_loo)
  gcv_obj <- lgspline:::.compute_gcvu(par, c(), env_gcv)
  gcv_grad <- lgspline:::.compute_gcvu_gradient(par, c(), env = env_gcv)

  expect_true(is.finite(loo_obj$criterion_value))
  expect_true(all(is.finite(loo_grad$gradient)))
  expect_true(is.finite(gcv_obj$criterion_value))
  expect_true(all(is.finite(gcv_grad$gradient)))
})

test_that("use_custom_bfgs = FALSE works for both tuning criteria", {
  t <- seq(-1, 1, length.out = 12)
  y <- t^2 + 0.1 * t

  fit_loo <- lgspline(
    cbind(t), y,
    K = 0,
    opt = TRUE,
    use_custom_bfgs = FALSE,
    tuning_criterion = "loo",
    initial_wiggle = c(1e-4, 5e-4),
    initial_flat = c(0.2, 0.5),
    unique_penalty_per_predictor = FALSE,
    unique_penalty_per_partition = FALSE,
    just_linear_without_interactions = 1,
    standardize_response = FALSE,
    include_warnings = FALSE
  )
  fit_gcv <- lgspline(
    cbind(t), y,
    K = 0,
    opt = TRUE,
    use_custom_bfgs = FALSE,
    tuning_criterion = "gcv",
    gcv_gamma = 1.4,
    initial_wiggle = c(1e-4, 5e-4),
    initial_flat = c(0.2, 0.5),
    unique_penalty_per_predictor = FALSE,
    unique_penalty_per_partition = FALSE,
    just_linear_without_interactions = 1,
    standardize_response = FALSE,
    include_warnings = FALSE
  )

  expect_true(all(is.finite(fit_loo$penalties$Lambda)))
  expect_true(all(is.finite(fit_gcv$penalties$Lambda)))
})

test_that("post-optimization sample-size adjustment shrinks tuned penalties for both criteria", {
  t <- seq(-1, 1, length.out = 14)
  y <- sin(pi * t) + 0.05 * t^2
  initial_wiggle <- c(1e-4, 5e-4)
  initial_flat <- c(0.2, 0.5)
  sample_size_adjustment <- (length(t) - 1) / (length(t) + 1)

  base_fit <- lgspline(
    cbind(t), y,
    K = 1,
    opt = FALSE,
    wiggle_penalty = 5e-3,
    flat_ridge_penalty = 0.3,
    unique_penalty_per_predictor = FALSE,
    unique_penalty_per_partition = FALSE,
    just_linear_without_interactions = 1,
    standardize_response = FALSE,
    include_warnings = FALSE
  )

  for (criterion in c("loo", "gcv")) {
    env <- build_tuning_env_from_fit(
      base_fit,
      tuning_criterion = criterion,
      gcv_gamma = 1.4
    )
    criterion_fxn <- if (criterion == "loo") {
      lgspline:::.compute_loocv
    } else {
      lgspline:::.compute_gcvu
    }
    best_start <- lgspline:::.tune_grid_search(
      log(initial_wiggle),
      log(initial_flat),
      c(),
      criterion_fxn,
      env,
      include_warnings = FALSE
    )
    opt_res <- optim(
      c(best_start, c()),
      fn = function(par) criterion_fxn(par, c(), env)$criterion_value,
      method = "BFGS"
    )
    raw_penalties <- exp(opt_res$par[1:2])

    fit_tuned <- lgspline(
      cbind(t), y,
      K = 1,
      opt = TRUE,
      use_custom_bfgs = FALSE,
      tuning_criterion = criterion,
      gcv_gamma = 1.4,
      initial_wiggle = initial_wiggle,
      initial_flat = initial_flat,
      unique_penalty_per_predictor = FALSE,
      unique_penalty_per_partition = FALSE,
      just_linear_without_interactions = 1,
      standardize_response = FALSE,
      include_warnings = FALSE
    )

    expect_equal(
      fit_tuned$penalties$wiggle_penalty,
      raw_penalties[1] * sample_size_adjustment,
      tolerance = 1e-8
    )
    expect_equal(
      fit_tuned$penalties$flat_ridge_penalty,
      raw_penalties[2] * sample_size_adjustment,
      tolerance = 1e-8
    )
    expect_lt(fit_tuned$penalties$wiggle_penalty, raw_penalties[1])
    expect_lt(fit_tuned$penalties$flat_ridge_penalty, raw_penalties[2])
  }
})

test_that(".damped_bfgs never returns a GCV solution worse than the grid start", {
  set.seed(1234)
  t <- seq(-2, 2, length.out = 30)
  y <- sin(2 * t) + 0.15 * t^2

  fit <- lgspline(
    cbind(t), y,
    K = 1,
    opt = FALSE,
    wiggle_penalty = 5e-3,
    flat_ridge_penalty = 0.3,
    unique_penalty_per_predictor = FALSE,
    unique_penalty_per_partition = FALSE,
    just_linear_without_interactions = 1,
    standardize_response = FALSE,
    include_warnings = FALSE
  )

  env <- build_tuning_env_from_fit(fit, tuning_criterion = "gcv", gcv_gamma = 1.4)
  log_penalty_vec <- c()
  log_initial_wiggle <- log(c(1e-4, 5e-4, 5e-3))
  log_initial_flat <- log(c(0.2, 0.5, 1))

  best_start <- lgspline:::.tune_grid_search(
    log_initial_wiggle,
    log_initial_flat,
    log_penalty_vec,
    lgspline:::.compute_gcvu,
    env,
    include_warnings = FALSE
  )

  start_value <- lgspline:::.compute_gcvu(best_start, log_penalty_vec, env)$criterion_value
  res <- lgspline:::.damped_bfgs(
    c(best_start, log_penalty_vec),
    log_penalty_vec,
    lgspline:::.compute_gcvu,
    lgspline:::.compute_gcvu_gradient,
    env,
    tol = 1e-6
  )

  expect_lte(res$criterion_value, start_value + 1e-10)
})

test_that(".damped_bfgs preserves the start point when the first step is worse", {
  env <- list()
  start <- c(0, 0)

  criterion_fxn <- function(par, log_penalty_vec, env, ...) {
    list(criterion_value = sum(par^2))
  }
  gr_fxn <- function(par, log_penalty_vec, outlist, env, ...) {
    list(
      gradient = c(-10, -10),
      outlist = list(criterion_value = sum(par^2))
    )
  }

  res <- lgspline:::.damped_bfgs(
    par = start,
    log_penalty_vec = c(),
    criterion_fxn = criterion_fxn,
    gr_fxn = gr_fxn,
    env = env,
    tol = 1e-8,
    max_iter = 2
  )

  expect_equal(res$par, start, tolerance = 1e-12)
  expect_equal(res$criterion_value, 0, tolerance = 1e-12)
})

test_that(".expand_tuning_grid adds extra candidates only once workers exceed six", {
  initial_grid <- expand.grid(wiggle = log(c(1e-4, 1e-2)),
                              flat = log(c(0.5, 5)))

  expanded <- lgspline:::.expand_tuning_grid(
    initial_grid,
    log_initial_wiggle = log(c(1e-4, 1e-2)),
    log_initial_flat = log(c(0.5, 5)),
    n_workers = 7L,
    parallel_grideval = TRUE
  )

  unchanged <- lgspline:::.expand_tuning_grid(
    initial_grid,
    log_initial_wiggle = log(c(1e-4, 1e-2)),
    log_initial_flat = log(c(0.5, 5)),
    n_workers = 6L,
    parallel_grideval = TRUE
  )

  expect_equal(unchanged, initial_grid)
  expect_equal(expanded[seq_len(nrow(initial_grid)), ], initial_grid)
  expect_equal(nrow(expanded), nrow(initial_grid) + 1L)
  expect_true(all(exp(expanded$wiggle[-seq_len(nrow(initial_grid))]) >= 1e-5))
  expect_true(all(exp(expanded$wiggle[-seq_len(nrow(initial_grid))]) <= 1e-1))
  expect_true(all(exp(expanded$flat[-seq_len(nrow(initial_grid))]) >= 0.05))
  expect_true(all(exp(expanded$flat[-seq_len(nrow(initial_grid))]) <= 50))
})

test_that("parallel grid evaluation matches sequential grid evaluation when no expansion is used", {
  case <- make_gaussian_linear_fit(wiggle_penalty = 5e-3, flat_ridge_penalty = 0.3)
  env <- build_tuning_env_from_fit(case$fit, tuning_criterion = "gcv", gcv_gamma = 1.4)
  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  lib_dir <- normalizePath(.libPaths()[1], winslash = "/", mustWork = TRUE)
  parallel::clusterExport(cl, "lib_dir", envir = environment())
  invisible(parallel::clusterEvalQ(cl, {
    .libPaths(c(lib_dir, .libPaths()))
    suppressPackageStartupMessages(library(lgspline))
    NULL
  }))

  seq_best <- lgspline:::.tune_grid_search(
    log(c(1e-4, 5e-4, 5e-3)),
    log(c(0.2, 0.5, 1)),
    c(),
    lgspline:::.compute_gcvu,
    env,
    include_warnings = FALSE,
    parallel_grideval = FALSE,
    cl = NULL
  )

  par_best <- lgspline:::.tune_grid_search(
    log(c(1e-4, 5e-4, 5e-3)),
    log(c(0.2, 0.5, 1)),
    c(),
    lgspline:::.compute_gcvu,
    env,
    include_warnings = FALSE,
    parallel_grideval = TRUE,
    cl = cl
  )

  expect_equal(par_best, seq_best, tolerance = 1e-10)
})

test_that("parallel damped BFGS returns a finite tuning solution", {
  set.seed(1234)
  t <- seq(-2, 2, length.out = 30)
  y <- sin(2 * t) + 0.15 * t^2

  fit <- lgspline(
    cbind(t), y,
    K = 1,
    opt = FALSE,
    wiggle_penalty = 5e-3,
    flat_ridge_penalty = 0.3,
    unique_penalty_per_predictor = FALSE,
    unique_penalty_per_partition = FALSE,
    just_linear_without_interactions = 1,
    standardize_response = FALSE,
    include_warnings = FALSE
  )

  env <- build_tuning_env_from_fit(fit, tuning_criterion = "gcv", gcv_gamma = 1.4)
  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  lib_dir <- normalizePath(.libPaths()[1], winslash = "/", mustWork = TRUE)
  parallel::clusterExport(cl, "lib_dir", envir = environment())
  invisible(parallel::clusterEvalQ(cl, {
    .libPaths(c(lib_dir, .libPaths()))
    suppressPackageStartupMessages(library(lgspline))
    NULL
  }))
  env$cl <- cl

  log_penalty_vec <- c()
  best_start <- lgspline:::.tune_grid_search(
    log(c(1e-4, 5e-4, 5e-3)),
    log(c(0.2, 0.5, 1)),
    log_penalty_vec,
    lgspline:::.compute_gcvu,
    env,
    include_warnings = FALSE,
    parallel_grideval = TRUE,
    cl = cl
  )

  start_value <- lgspline:::.compute_gcvu(best_start, log_penalty_vec, env)$criterion_value
  res <- lgspline:::.damped_bfgs(
    c(best_start, log_penalty_vec),
    log_penalty_vec,
    lgspline:::.compute_gcvu,
    lgspline:::.compute_gcvu_gradient,
    env,
    tol = 1e-6,
    parallel_bfgs = TRUE
  )

  expect_true(all(is.finite(res$par)))
  expect_true(is.finite(res$criterion_value))
  expect_lte(res$criterion_value, start_value + 1e-10)
})

test_that("adaptive tuning tolerance remains small and scale-aware", {
  tol <- 10 * sqrt(.Machine$double.eps)

  expect_equal(lgspline:::.tuning_adaptive_criterion_tol(tol, 1e-4),
               tol)
  expect_equal(lgspline:::.tuning_adaptive_criterion_tol(tol, 0.27),
               2.7e-6,
               tolerance = 1e-12)
  expect_equal(lgspline:::.tuning_adaptive_criterion_tol(tol, 10),
               1e-5,
               tolerance = 1e-12)
})

test_that("tuning active-set warm start keeps valid inequality columns", {
  env <- list(
    qp_Amat = matrix(0, nrow = 3, ncol = 5),
    active_set_cache = new.env(parent = emptyenv())
  )

  lgspline:::.tuning_store_active_set(
    env, list(active_ineq = c(2, 4, 9, NA, 2))
  )
  expect_equal(lgspline:::.tuning_seed_active_set(env), c(2L, 4L))

  lgspline:::.tuning_clear_active_set(env)
  expect_equal(lgspline:::.tuning_seed_active_set(env), integer(0))
})

test_that("custom GCV tuning keeps partition penalties in range", {
  set.seed(1234)
  t <- sort(runif(80, -3, 3))
  eta <- sin(t) - 0.2 * t^2
  y <- rbinom(length(t), 1, plogis(eta))
  dat <- data.frame(t = t, y = y)

  fit_fd <- lgspline(
    y ~ spl(t),
    K = 1,
    dat,
    tuning_criterion = "gcv",
    use_custom_bfgs = FALSE,
    unique_penalty_per_predictor = FALSE,
    unique_penalty_per_partition = TRUE,
    include_warnings = FALSE
  )
  fit_custom <- lgspline(
    y ~ spl(t),
    K = 1,
    dat,
    tuning_criterion = "gcv",
    use_custom_bfgs = TRUE,
    unique_penalty_per_predictor = FALSE,
    unique_penalty_per_partition = TRUE,
    include_warnings = FALSE
  )

  ratio <- max(fit_custom$penalties$other_penalties) /
    max(fit_fd$penalties$other_penalties)

  expect_lte(ratio, 10)
})
