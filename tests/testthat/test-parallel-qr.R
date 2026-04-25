expect_flag <- function(flag, message) {
  if (!isTRUE(flag)) testthat::fail(message)
}


test_that("parallel_qr least-squares helper matches lm.fit", {
  skip_on_cran()

  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl))

  set.seed(1234)
  X <- cbind(1, matrix(rnorm(600 * 11), nrow = 600))
  y <- rnorm(600)

  fit_serial <- .lm.fit(X, y)
  fit_parallel <- lgspline:::.parallel_qr_lm_fit(
    X = X,
    y = y,
    parallel_qr = TRUE,
    cl = cl
  )

  expect_equal(c(fit_parallel$coefficients),
               c(fit_serial$coefficients),
               tolerance = 1e-7)
  expect_equal(c(fit_parallel$residuals),
               c(fit_serial$residuals),
               tolerance = 1e-7)
})


test_that("parallel_qr constraint reduction preserves the column space", {
  skip_on_cran()

  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl))

  set.seed(4321)
  A_base <- matrix(rnorm(300 * 6), nrow = 300)
  A <- cbind(A_base,
             A_base[, 1] + 2 * A_base[, 2],
             A_base[, 3])

  A_serial <- lgspline:::.reduce_constraint_basis(A)
  A_parallel <- lgspline:::.reduce_constraint_basis(
    A = A,
    parallel_qr = TRUE,
    cl = cl
  )

  P_serial <- A_serial %*% t(A_serial)
  P_parallel <- A_parallel %*% t(A_parallel)

  expect_equal(P_parallel, P_serial, tolerance = 1e-6)
})


test_that("parallel_qr projection helper matches the serial path", {
  skip_on_cran()

  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl))

  set.seed(1234)
  K <- 80
  p_expansions <- 4
  P <- p_expansions * (K + 1)
  A <- matrix(rnorm(P * 12), nrow = P)
  Ghalf <- lapply(1:(K + 1), function(k) {
    diag(runif(p_expansions, 0.75, 1.25))
  })
  GhalfXy <- cbind(rnorm(P))

  fit_serial <- lgspline:::.lagrangian_project(
    GhalfXy = GhalfXy,
    Ghalf = Ghalf,
    A = A,
    K = K,
    p_expansions = p_expansions,
    R_constraints = ncol(A),
    constraint_value_vectors = list(),
    family = gaussian(),
    parallel_aga = FALSE,
    parallel_matmult = FALSE,
    cl = NULL,
    chunk_size = 1,
    num_chunks = 0,
    rem_chunks = 0,
    parallel_qr = FALSE
  )
  fit_parallel <- lgspline:::.lagrangian_project(
    GhalfXy = GhalfXy,
    Ghalf = Ghalf,
    A = A,
    K = K,
    p_expansions = p_expansions,
    R_constraints = ncol(A),
    constraint_value_vectors = list(),
    family = gaussian(),
    parallel_aga = FALSE,
    parallel_matmult = FALSE,
    cl = cl,
    chunk_size = 1,
    num_chunks = 0,
    rem_chunks = 0,
    parallel_qr = TRUE
  )

  expect_equal(unlist(fit_parallel),
               unlist(fit_serial),
               tolerance = 1e-6)
})


test_that("QP QR pivot keeps equality and cross-partition columns", {
  qp_Amat <- cbind(
    c(1, 0, 0, 0),  # equality
    c(1, 2, 0, 0),  # partition 1
    c(2, 4, 0, 0),  # dependent partition 1
    c(0, 0, 1, 0),  # partition 2
    c(0, 1, 1, 0)   # cross-partition, keep as-is
  )
  qp_bvec <- c(0, 0, 0, 0, 0)

  out <- lgspline:::.reduce_qp_constraint_block(
    qp_Amat = qp_Amat,
    qp_bvec = qp_bvec,
    qp_meq = 1,
    p_expansions = 2,
    K = 1
  )

  expect_equal(out$meq, 1L)
  expect_equal(ncol(out$Amat), 4L)
  expect_equal(out$Amat[, 1], qp_Amat[, 1])
  expect_flag(any(apply(out$Amat, 2, function(col) identical(col, qp_Amat[, 4]))),
              "Partition-2 local QP column was dropped unexpectedly.")
  expect_flag(any(apply(out$Amat, 2, function(col) identical(col, qp_Amat[, 5]))),
              "Cross-partition QP column was dropped unexpectedly.")
})


test_that("process_qp can thin derivative constraints with QP QR pivoting", {
  set.seed(20260424)

  t <- seq(-3, 3, length.out = 120)
  y <- exp(0.25 * t) + rnorm(length(t), 0, 0.02)
  fit <- lgspline::lgspline(
    cbind(t = t),
    y,
    K = 2,
    opt = FALSE,
    standardize_response = FALSE,
    qp_positive_derivative = "t",
    include_cubic_terms = FALSE,
    include_quartic_terms = FALSE,
    include_constrain_second_deriv = FALSE,
    include_warnings = FALSE
  )

  fit_obj <- if (!is.null(fit$model_fit)) fit$model_fit else fit
  qp_full <- lgspline:::process_qp(
    X = fit_obj$X,
    K = fit_obj$K,
    p_expansions = fit_obj$p,
    order_list = fit_obj$order_list,
    colnm_expansions = fit_obj$raw_expansion_names,
    expansion_scales = fit_obj$expansion_scales,
    power1_cols = fit_obj$power1_cols,
    power2_cols = fit_obj$power2_cols,
    nonspline_cols = fit_obj$nonspline_cols,
    interaction_single_cols = fit_obj$interaction_single_cols,
    interaction_quad_cols = fit_obj$interaction_quad_cols,
    triplet_cols = fit_obj$triplet_cols,
    include_2way_interactions =
      fit_obj$.fit_call_args$include_2way_interactions,
    include_3way_interactions =
      fit_obj$.fit_call_args$include_3way_interactions,
    include_quadratic_interactions =
      fit_obj$.fit_call_args$include_quadratic_interactions,
    family = fit_obj$family,
    mean_y = fit_obj$mean_y,
    sd_y = fit_obj$sd_y,
    N_obs = fit_obj$N,
    qp_positive_derivative = 1,
    og_cols = fit_obj$og_cols,
    qr_pivot_inequality_constraints = FALSE,
    include_warnings = FALSE
  )
  qp_pivot <- lgspline:::process_qp(
    X = fit_obj$X,
    K = fit_obj$K,
    p_expansions = fit_obj$p,
    order_list = fit_obj$order_list,
    colnm_expansions = fit_obj$raw_expansion_names,
    expansion_scales = fit_obj$expansion_scales,
    power1_cols = fit_obj$power1_cols,
    power2_cols = fit_obj$power2_cols,
    nonspline_cols = fit_obj$nonspline_cols,
    interaction_single_cols = fit_obj$interaction_single_cols,
    interaction_quad_cols = fit_obj$interaction_quad_cols,
    triplet_cols = fit_obj$triplet_cols,
    include_2way_interactions =
      fit_obj$.fit_call_args$include_2way_interactions,
    include_3way_interactions =
      fit_obj$.fit_call_args$include_3way_interactions,
    include_quadratic_interactions =
      fit_obj$.fit_call_args$include_quadratic_interactions,
    family = fit_obj$family,
    mean_y = fit_obj$mean_y,
    sd_y = fit_obj$sd_y,
    N_obs = fit_obj$N,
    qp_positive_derivative = 1,
    og_cols = fit_obj$og_cols,
    qr_pivot_inequality_constraints = TRUE,
    include_warnings = FALSE
  )

  expect_flag(qp_full$quadprog,
              "Unpivoted derivative QP block was not marked active.")
  expect_flag(qp_pivot$quadprog,
              "Pivoted derivative QP block was not marked active.")
  expect_equal(qp_full$qp_meq, qp_pivot$qp_meq)
  expect_flag(ncol(qp_pivot$qp_Amat) <= ncol(qp_full$qp_Amat),
              "QP QR pivoting increased the number of derivative constraints.")
})
