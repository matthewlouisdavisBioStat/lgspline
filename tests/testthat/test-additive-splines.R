additive_term_coef_has <- function(fit, term, pattern) {
  nm <- names(unlist(fit$B))
  any(grepl(paste0("^smooth", term, "\\."), nm) & grepl(pattern, nm))
}


additive_coef_has <- function(fit, pattern) {
  any(grepl(pattern, names(unlist(fit$B))))
}


additive_term_coef_absmax <- function(fit, term, pattern) {
  vals <- unlist(fit$B)
  ind <- grepl(paste0("^smooth", term, "\\."), names(vals)) &
    grepl(pattern, names(vals))
  if(!any(ind)) return(0)
  max(abs(vals[ind]))
}


additive_zero_anchor_cols <- function(fit, term) {
  term_fit <- fit$additive_terms[[term]]
  cols <- term_fit$additive_zero_main_cols
  if(length(cols) == 0L) return(TRUE)
  vals <- unlist(term_fit$B, use.names = FALSE)
  p <- term_fit$p
  ind <- unlist(lapply(seq_len(term_fit$K + 1L), function(k) {
    (k - 1L) * p + cols
  }), use.names = FALSE)
  max(abs(vals[ind])) < 1e-6
}


additive_logistic_link_did <- function(fit) {
  nd <- data.frame(t1 = c(1, -1, 1, -1),
                   t2 = 0,
                   z = c(1, 1, 0, 0))
  mu <- c(predict(fit, nd))
  mu <- pmin(pmax(mu, 1e-8), 1 - 1e-8)
  sum(c(1, -1, -1, 1) * fit$family$linkfun(mu))
}


additive_equalities_touch_offset <- function(fit) {
  eq_fun <- getFromNamespace(".additive_combined_equalities", "lgspline")
  eq <- eq_fun(fit$additive_terms)
  if(ncol(eq$Amat) == 0L) return(FALSE)
  rows <- grepl(".additive_offset", names(unlist(fit$B)), fixed = TRUE)
  any(colSums(abs(eq$Amat[rows, , drop = FALSE])) > 1e-10)
}


additive_system_offset_design_zero <- function(fit) {
  system_fun <- getFromNamespace(".additive_combined_system", "lgspline")
  system <- system_fun(fit)
  cols <- grepl(".additive_offset", names(unlist(fit$B)), fixed = TRUE)
  if(!any(cols)) return(TRUE)
  max(abs(system$X[, cols, drop = FALSE])) < 1e-10
}


test_that("separate spl() calls default to additive spline groups", {
  set.seed(101)
  n <- 50
  dat <- data.frame(
    t1 = runif(n, -2, 2),
    t2 = runif(n, -1, 1)
  )
  dat$y <- sin(dat$t1) + dat$t2^2 + rnorm(n, sd = 0.05)

  fit_add <- lgspline(
    y ~ spl(t1) + spl(t2),
    dat,
    K = 1,
    opt = FALSE,
    return_varcovmat = FALSE,
    include_warnings = FALSE,
    use_s_alias = TRUE,
    additive_max_iter = 3
  )
  fit_s_alias <- lgspline(
    y ~ s(t1) + s(t2),
    dat,
    K = 1,
    opt = FALSE,
    return_varcovmat = FALSE,
    include_warnings = FALSE,
    additive_max_iter = 3
  )
  fit_joined <- lgspline(
    y ~ spl(t1, t2),
    dat,
    K = 1,
    opt = FALSE,
    return_varcovmat = FALSE,
    include_warnings = FALSE
  )

  expect_s3_class(fit_add, "additive_lgspline")
  expect_s3_class(fit_s_alias, "additive_lgspline")
  expect_s3_class(fit_joined, "lgspline")
  expect_false(inherits(fit_joined, "additive_lgspline"))
  expect_equal(length(fit_add$additive_terms), 2)
  expect_equal(length(fit_s_alias$additive_terms), 2)
  expect_equal(length(fit_add$spline_groups), 2)
  expect_equal(fit_s_alias$spline_groups, fit_add$spline_groups)
  expect_equal(unname(fit_add$K), c(1, 1))
  expect_equal(length(predict(fit_add)), n)
  expect_equal(c(predict(fit_s_alias)), c(predict(fit_add)), tolerance = 1e-8)
})


test_that("additive logistic fits preserve anchor interaction signal", {
  set.seed(115)
  n <- 500
  dat <- data.frame(
    t1 = runif(n, -1.5, 1.5),
    t2 = runif(n, -1.5, 1.5),
    z = rbinom(n, 1, 0.5)
  )
  eta <- -0.60 + 0.55 * sin(pi * dat$t1) -
    0.45 * cos(pi * dat$t2 / 2) + 0.35 * dat$z +
    0.9 * dat$t1 * dat$z
  dat$y <- rbinom(n, 1, plogis(eta))

  fit <- lgspline(
    y ~ spl(t1) + spl(t2) + z + t1:z,
    dat,
    family = binomial(),
    K = 1,
    opt = FALSE,
    return_varcovmat = FALSE,
    include_warnings = FALSE,
    additive_max_iter = 4
  )

  expect_s3_class(fit, "additive_lgspline")
  expect_false(additive_equalities_touch_offset(fit))
  expect_true(additive_system_offset_design_zero(fit))
  expect_gt(additive_logistic_link_did(fit), 0.75)
  expect_true(additive_term_coef_absmax(fit, 1, "t1xz") > 0.35)
})


test_that("additive formula interactions are stable across term permutations", {
  set.seed(112)
  n <- 48
  dat <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1),
    x3 = runif(n, -1, 1),
    z = rbinom(n, 1, 0.45),
    w = rnorm(n),
    f = factor(rep(letters[1:3], length.out = n))
  )
  dat$off <- rnorm(n, sd = 0.05)
  dat$y <- sin(dat$x1) + dat$x2^2 - 0.25 * dat$x3 +
    0.4 * dat$x1 * dat$z - 0.3 * dat$x2 * dat$w +
    c(a = -0.2, b = 0.1, c = 0.25)[dat$f] + dat$off +
    rnorm(n, sd = 0.06)

  fit0 <- function(formula) {
    lgspline(
      formula,
      dat,
      K = 1,
      opt = FALSE,
      return_varcovmat = FALSE,
      include_warnings = FALSE,
      additive_max_iter = 2
    )
  }

  cases <- list(
    list(y ~ spl(x1) + spl(x2) + x1:z, 1, "x1xz", 2),
    list(y ~ spl(x2) + spl(x1) + x1:z, 2, "x1xz", 2),
    list(y ~ spl(x1) + spl(x2) + x2:z, 2, "x2xz", 2),
    list(y ~ spl(x2) + spl(x1) + x2:z, 1, "x2xz", 2),
    list(y ~ spl(x2) + spl(x1) + z * x1, 2, "x1xz", 2),
    list(y ~ s(x1) + s(x2) + z:x2, 2, "x2xz", 2),
    list(y ~ 0 + spl(x2) + spl(x1) + x1:z + offset(off), 2,
         "x1xz", 2),
    list(y ~ spl(x2) + spl(x1) + x1:f, 2, "x1x", NA_integer_),
    list(y ~ spl(x3) + spl(x2) + spl(x1) + x1:z + x2:w, 3,
         "x1xz", 2)
  )

  for(case in cases) {
    fit <- fit0(case[[1]])
    expect_s3_class(fit, "additive_lgspline")
    expect_true(additive_term_coef_has(fit, case[[2]], case[[3]]))
    expect_true(additive_coef_has(fit, case[[3]]))
    if(!is.na(case[[4]])) {
      expect_true(additive_zero_anchor_cols(fit, case[[4]]))
    }
    expect_true(all(is.finite(predict(fit))))
  }
})


test_that("direct additive fits keep anchor-only interactions identified", {
  set.seed(113)
  n <- 46
  predictors <- cbind(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1),
    z = rbinom(n, 1, 0.5),
    w = rnorm(n)
  )
  y <- sin(predictors[, "x1"]) + predictors[, "x2"]^2 +
    0.5 * predictors[, "x1"] * predictors[, "z"] +
    rnorm(n, sd = 0.06)

  fit <- lgspline.fit(
    predictors,
    y,
    K = 1,
    opt = FALSE,
    return_varcovmat = FALSE,
    spline_groups = list("x2", "x1"),
    just_linear_with_interactions = c("z", "w"),
    include_warnings = FALSE,
    additive_max_iter = 2
  )
  fit_excl <- lgspline.fit(
    predictors,
    y,
    K = 1,
    opt = FALSE,
    return_varcovmat = FALSE,
    spline_groups = list("x2", "x1"),
    just_linear_with_interactions = c("z", "w"),
    exclude_these_expansions = "_3_x_4_",
    include_warnings = FALSE,
    additive_max_iter = 2
  )

  expect_s3_class(fit, "additive_lgspline")
  expect_true(additive_term_coef_has(fit, 2, "x1xz"))
  expect_true(additive_term_coef_has(fit, 1, "x2xw"))
  expect_true(additive_term_coef_has(fit, 1, "\\.zxw$"))
  expect_lt(additive_term_coef_absmax(fit, 2, "\\.zxw$"), 1e-6)
  expect_false(additive_coef_has(fit_excl, "\\.zxw$"))
})


test_that("additive interactions work with GLM, blockfit, QP, and correlation", {
  set.seed(114)
  n <- 42
  dat <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1),
    z = rbinom(n, 1, 0.45),
    f = factor(rep(letters[1:3], length.out = n))
  )
  eta <- -0.15 + 0.55 * dat$x1 - 0.25 * dat$x2 +
    0.35 * dat$x1 * dat$z +
    c(a = -0.1, b = 0.1, c = 0.2)[dat$f]
  dat$y <- rbinom(n, 1, plogis(eta))
  V <- 0.25 ^ abs(outer(seq_len(n), seq_len(n), "-"))

  fit <- lgspline(
    y ~ spl(x2) + spl(x1) + x1:z + f,
    dat,
    family = binomial(),
    K = 1,
    opt = FALSE,
    blockfit = TRUE,
    return_varcovmat = FALSE,
    VhalfInv = matinvsqrt(V),
    qp_positive_derivative = "x1",
    include_warnings = FALSE,
    additive_max_iter = 2
  )

  expect_s3_class(fit, "additive_lgspline")
  expect_true(additive_term_coef_has(fit, 2, "x1xz"))
  expect_true(all(is.finite(predict(fit))))
  expect_false(identical(fit$additive_terms[[2]]$quadprog_list, list(NA)))
  expect_true(length(fit$additive_terms[[1]]$constraint_vectors) > 0)
})


test_that("direct lgspline.fit accepts additive spline_groups", {
  set.seed(102)
  n <- 45
  predictors <- cbind(
    t1 = runif(n, -2, 2),
    t2 = runif(n, -1, 1)
  )
  y <- cos(predictors[, "t1"]) + predictors[, "t2"]^2 + rnorm(n, sd = 0.05)

  fit <- lgspline.fit(
    predictors,
    y,
    K = 1,
    opt = FALSE,
    return_varcovmat = FALSE,
    spline_groups = list("t1", "t2"),
    include_warnings = FALSE,
    additive_max_iter = 3
  )

  expect_s3_class(fit, "additive_lgspline")
  expect_equal(length(fit$additive_terms), 2)
  expect_equal(length(predict(fit, predictors[1:4, ])), 4)
})


test_that("additive splines keep predictor-specific QP constraints", {
  set.seed(103)
  n <- 48
  dat <- data.frame(
    t1 = sort(runif(n, 0, 2)),
    t2 = runif(n, -1, 1)
  )
  dat$y <- dat$t1 + dat$t2^2 + rnorm(n, sd = 0.04)

  fit <- lgspline(
    y ~ spl(t1) + spl(t2),
    dat,
    K = 1,
    opt = FALSE,
    return_varcovmat = FALSE,
    qp_positive_derivative = "t1",
    qp_observations = list("t1:qp_positive_derivative" = 1:20),
    include_warnings = FALSE,
    additive_max_iter = 3
  )

  expect_s3_class(fit, "additive_lgspline")
  expect_false(is.null(fit$additive_terms[[1]]$quadprog_list$qp_Amat))
  expect_true(is.null(fit$additive_terms[[2]]$qp_info) ||
                is.null(fit$additive_terms[[2]]$quadprog_list$qp_Amat) ||
                identical(fit$additive_terms[[2]]$quadprog_list, list(NA)))
})


test_that("additive splines accept custom and per-term QP constraints", {
  set.seed(104)
  n <- 40
  predictors <- cbind(
    t1 = runif(n, -1, 1),
    t2 = runif(n, -1, 1)
  )
  y <- predictors[, "t1"]^2 - predictors[, "t2"] + rnorm(n, sd = 0.05)

  qp_Amat_fxn <- function(N, p, K, X, colnm, scales, deriv_fxn, ...) {
    out <- matrix(0, p * (K + 1), 1)
    out[1, 1] <- 1
    out
  }
  qp_bvec_fxn <- function(qp_Amat, ...) rep(-1e6, ncol(qp_Amat))
  qp_meq_fxn <- function(qp_Amat, ...) 0

  fit_fxn <- lgspline.fit(
    predictors,
    y,
    K = 1,
    opt = FALSE,
    return_varcovmat = FALSE,
    spline_groups = list("t1", "t2"),
    qp_Amat_fxn = qp_Amat_fxn,
    qp_bvec_fxn = qp_bvec_fxn,
    qp_meq_fxn = qp_meq_fxn,
    include_warnings = FALSE,
    additive_max_iter = 2
  )

  dummy_term <- lgspline.fit(
    cbind(t = predictors[, 1], .additive_offset = 0),
    y,
    K = 1,
    opt = FALSE,
    dummy_fit = TRUE,
    return_varcovmat = FALSE,
    offset = 2,
    just_linear_without_interactions = 2,
    do_not_cluster_on_these = 2,
    include_warnings = FALSE
  )
  term_qp <- matrix(0, dummy_term$P, 1)
  term_qp[1, 1] <- 1

  fit_list <- lgspline.fit(
    predictors,
    y,
    K = 1,
    opt = FALSE,
    return_varcovmat = FALSE,
    spline_groups = list("t1", "t2"),
    qp_Amat = list(term_qp, term_qp),
    qp_bvec = list(-1e6, -1e6),
    qp_meq = list(0, 0),
    include_warnings = FALSE,
    additive_max_iter = 2
  )

  expect_s3_class(fit_fxn, "additive_lgspline")
  expect_s3_class(fit_list, "additive_lgspline")
  expect_false(identical(fit_fxn$additive_terms[[1]]$quadprog_list, list(NA)))
  expect_false(identical(fit_list$additive_terms[[1]]$quadprog_list, list(NA)))
  expect_error(
    lgspline.fit(
      predictors,
      y,
      K = 1,
      opt = FALSE,
      spline_groups = list("t1", "t2"),
      qp_Amat = term_qp,
      qp_bvec = -1e6,
      qp_meq = 0,
      include_warnings = FALSE
    ),
    "Global pre-built"
  )
})


test_that("additive splines preserve factor constraints under blockfit", {
  set.seed(105)
  n <- 45
  dat <- data.frame(
    t1 = runif(n, -1, 1),
    t2 = runif(n, -1, 1),
    grp = factor(rep(letters[1:3], length.out = n))
  )
  grp_eff <- c(a = -0.2, b = 0.1, c = 0.3)
  dat$y <- sin(dat$t1) + dat$t2^2 + grp_eff[dat$grp] +
    rnorm(n, sd = 0.05)

  fit <- lgspline(
    y ~ spl(t1) + spl(t2) + grp,
    dat,
    K = 1,
    opt = FALSE,
    blockfit = TRUE,
    return_varcovmat = FALSE,
    include_warnings = FALSE,
    additive_max_iter = 2
  )

  expect_s3_class(fit, "additive_lgspline")
  expect_true(any(vapply(fit$additive_terms, function(term) {
    isTRUE(term$blockfit_used)
  }, logical(1))))
  expect_true(length(fit$additive_terms[[1]]$constraint_vectors) > 0)
})


test_that("additive splines work with family objects and fixed correlation", {
  set.seed(106)
  n <- 42
  dat <- data.frame(
    t1 = runif(n, -1, 1),
    t2 = runif(n, -1, 1)
  )
  eta <- 0.2 + 0.5 * dat$t1 - 0.25 * dat$t2
  dat$y <- rpois(n, exp(eta))

  fit_glm <- lgspline(
    y ~ spl(t1) + spl(t2),
    dat,
    family = quasipoisson(),
    K = 1,
    opt = FALSE,
    return_varcovmat = FALSE,
    include_warnings = FALSE,
    additive_max_iter = 3
  )

  V <- diag(n)
  V[1:6, 1:6] <- 0.2
  diag(V) <- 1
  fit_cor <- lgspline(
    y ~ spl(t1) + spl(t2),
    dat,
    K = 1,
    opt = FALSE,
    return_varcovmat = FALSE,
    VhalfInv = matinvsqrt(V),
    include_warnings = FALSE,
    additive_max_iter = 3
  )

  expect_s3_class(fit_glm, "additive_lgspline")
  expect_true(all(predict(fit_glm) > 0))
  expect_s3_class(fit_cor, "additive_lgspline")
  expect_false(is.null(fit_cor$VhalfInv))
})


test_that("additive splines handle sparse and dense fixed correlation paths", {
  set.seed(107)
  n <- 36
  dat <- data.frame(
    t1 = runif(n, -1, 1),
    t2 = runif(n, -1, 1)
  )
  dat$y <- sin(dat$t1) - dat$t2 + rnorm(n, sd = 0.05)

  id <- rep(seq_len(6), each = 6)
  V_sparse <- diag(n)
  for(cl in unique(id)){
    ind <- which(id == cl)
    V_sparse[ind, ind] <- 0.25
    diag(V_sparse)[ind] <- 1
  }
  idx <- seq_len(n)
  V_dense <- 0.35 ^ abs(outer(idx, idx, "-"))

  fit_sparse <- lgspline(
    y ~ spl(t1) + spl(t2),
    dat,
    K = 1,
    opt = FALSE,
    return_varcovmat = FALSE,
    VhalfInv = matinvsqrt(V_sparse),
    include_warnings = FALSE,
    additive_max_iter = 2
  )
  fit_dense <- lgspline(
    y ~ spl(t1) + spl(t2),
    dat,
    K = 1,
    opt = FALSE,
    return_varcovmat = FALSE,
    VhalfInv = matinvsqrt(V_dense),
    include_warnings = FALSE,
    additive_max_iter = 2
  )

  expect_s3_class(fit_sparse, "additive_lgspline")
  expect_s3_class(fit_dense, "additive_lgspline")
  expect_true(all(is.finite(predict(fit_sparse))))
  expect_true(all(is.finite(predict(fit_dense))))
})


test_that("additive splines handle estimated correlation structures", {
  set.seed(108)
  n <- 30
  dat <- data.frame(
    t1 = runif(n, -1, 1),
    t2 = runif(n, -1, 1),
    id = rep(seq_len(10), each = 3)
  )
  dat$y <- sin(dat$t1) + dat$t2 + rnorm(n, sd = 0.05)

  fit <- lgspline(
    y ~ spl(t1) + spl(t2),
    dat,
    K = 1,
    opt = FALSE,
    return_varcovmat = FALSE,
    correlation_id = dat$id,
    correlation_structure = "exchangeable",
    include_warnings = FALSE,
    additive_max_iter = 2
  )

  expect_s3_class(fit, "additive_lgspline")
  expect_equal(length(fit$VhalfInv_params_estimates), 1)
  expect_true(all(is.finite(predict(fit))))
  post <- generate_posterior(
    fit,
    draw_correlation = TRUE,
    num_draws = 2,
    include_warnings = FALSE
  )
  expect_equal(length(post$post_draw_correlation_params), 2)
  expect_equal(length(post$post_draw_coefficients), 2)
})


test_that("additive splines handle GLM, indicators, QP, and correlation together", {
  set.seed(109)
  n <- 36
  dat <- data.frame(
    t1 = runif(n, -1, 1),
    t2 = runif(n, -1, 1),
    grp = factor(rep(letters[1:3], length.out = n))
  )
  eta <- 0.2 + 0.4 * dat$t1 - 0.2 * dat$t2 +
    c(a = -0.1, b = 0.15, c = 0.25)[dat$grp]
  dat$y <- rpois(n, exp(eta))

  V <- 0.2 ^ abs(outer(seq_len(n), seq_len(n), "-"))
  fit <- lgspline(
    y ~ spl(t1) + spl(t2) + grp,
    dat,
    family = quasipoisson(),
    K = 1,
    opt = FALSE,
    blockfit = TRUE,
    return_varcovmat = FALSE,
    VhalfInv = matinvsqrt(V),
    qp_positive_derivative = "t1",
    include_warnings = FALSE,
    additive_max_iter = 2
  )

  expect_s3_class(fit, "additive_lgspline")
  expect_true(all(is.finite(predict(fit))))
  expect_false(identical(fit$additive_terms[[1]]$quadprog_list, list(NA)))
  expect_true(length(fit$additive_terms[[1]]$constraint_vectors) > 0)
})


test_that("additive post-fit methods work for three smooth terms", {
  set.seed(110)
  n <- 44
  dat <- data.frame(
    t1 = runif(n, -1, 1),
    t2 = runif(n, -1, 1),
    t3 = runif(n, -1, 1)
  )
  dat$y <- sin(dat$t1) + dat$t2^2 - 0.4 * dat$t3 +
    rnorm(n, sd = 0.05)

  fit <- lgspline(
    y ~ s(t1) + s(t2) + s(t3),
    dat,
    K = 1,
    opt = FALSE,
    return_varcovmat = TRUE,
    include_warnings = FALSE,
    additive_max_iter = 2
  )

  newdat <- dat[1:5, c("t1", "t2", "t3")]
  pred <- predict(fit, newdat)
  pred_se <- predict(fit, newdat, se.fit = TRUE)
  deriv <- predict(fit, newdat,
                  take_first_derivatives = TRUE,
                  take_second_derivatives = TRUE)
  printed <- capture.output(print(fit))
  sm <- summary(fit)
  printed_sm <- capture.output(print(sm))
  wald <- wald_univariate(fit)
  wald_sm <- capture.output(summary(wald))
  wald_coef <- coef(wald)
  wald_ci <- confint(wald)
  ci <- confint(fit)
  post <- generate_posterior(
    fit,
    new_predictors = newdat,
    include_posterior_predictive = TRUE,
    num_draws = 2,
    enforce_qp_constraints = TRUE
  )
  opt <- find_extremum(fit, vars = c("t1", "t2"), minimize = TRUE)
  eq <- NULL
  capture.output(
    eq <- equation(fit, first_derivative = "t1", show_bounds = FALSE)
  )
  loo <- leave_one_out(fit)
  integ <- integrate(fit, lower = -0.5, upper = 0.5,
                     vars = "t1", n_quad = 4)
  ll <- logLik(fit)
  pll <- prior_loglik(fit)

  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  plot(fit, vars = "t1", n_grid = 10)
  plot(wald, which = 1:3)
  grDevices::dev.off()

  expect_s3_class(fit, "additive_lgspline")
  expect_equal(length(fit$additive_terms), 3)
  expect_true(length(printed) > 0)
  expect_s3_class(sm, "summary.additive_lgspline")
  expect_true(length(printed_sm) > 0)
  expect_equal(length(pred), 5)
  expect_true(all(is.finite(pred)))
  expect_equal(length(pred_se$fit), 5)
  expect_true(all(is.finite(pred_se$se.fit)))
  expect_equal(length(deriv$first_deriv), 3)
  expect_true(is.matrix(wald$coefficients))
  expect_true(length(wald_sm) > 0)
  expect_true(length(wald_coef) > 0)
  expect_true(is.matrix(wald_ci))
  expect_true(is.matrix(ci))
  expect_equal(length(post$post_draw_coefficients), 2)
  expect_equal(nrow(post$post_pred_draw), 5)
  expect_equal(length(opt$t), 3)
  expect_true(is.list(eq))
  expect_equal(length(loo), n)
  expect_true(is.finite(integ))
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(pll))
})


test_that("additive post-fit methods work for four smooth terms", {
  set.seed(111)
  n <- 38
  dat <- data.frame(
    t1 = runif(n, -1, 1),
    t2 = runif(n, -1, 1),
    t3 = runif(n, -1, 1),
    t4 = runif(n, -1, 1)
  )
  dat$y <- 0.4 * dat$t1 - dat$t2^2 + sin(dat$t3) + 0.2 * dat$t4 +
    rnorm(n, sd = 0.06)

  fit <- lgspline(
    y ~ s(t1) + s(t2) + s(t3) + s(t4),
    dat,
    K = 1,
    opt = FALSE,
    return_varcovmat = TRUE,
    include_warnings = FALSE,
    additive_max_iter = 2
  )

  draw <- generate_posterior(fit, num_draws = 1)
  pred_draw <- predict(fit, dat[1:4, c("t1", "t2", "t3", "t4")],
                       B_predict = draw$post_draw_coefficients)
  ll_draw <- logLik(fit, B_predict = draw$post_draw_coefficients,
                    sigmasq_predict = draw$post_draw_sigmasq)
  no_corr <- try(
    generate_posterior_correlation(fit, include_warnings = FALSE),
    silent = TRUE
  )

  expect_s3_class(fit, "additive_lgspline")
  expect_equal(length(fit$additive_terms), 4)
  expect_equal(length(coef(fit)), 4)
  expect_equal(length(pred_draw), 4)
  expect_true(all(is.finite(pred_draw)))
  expect_s3_class(ll_draw, "logLik")
  expect_s3_class(no_corr, "try-error")
})
