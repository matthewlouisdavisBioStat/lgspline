test_that("Cox PH helpers", {

  ## Test 1: Partial log-likelihood at beta = 0 matches manual computation
  set.seed(1234)
  n <- 50
  time <- rexp(n, 1)
  status <- rbinom(n, 1, 0.7)
  ord <- order(time)
  time_s <- time[ord]
  status_s <- status[ord]

  event_times <- sort(unique(time_s[status_s == 1]))
  ll_manual <- 0
  for(tt in event_times) {
    Dg <- which(time_s == tt & status_s == 1)
    Rg <- which(time_s >= tt)
    ll_manual <- ll_manual + sum(0 - log(length(Rg)))
  }
  ll_fn <- loglik_cox(rep(0, n), status_s, time_s)
  expect_equal(ll_fn, ll_manual, tolerance = 1e-4)

  ## Test 2: Score at MLE should be approximately zero
  set.seed(1234)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  time <- rexp(n, exp(0.3 * x1 + 0.2 * x2))
  status <- rbinom(n, 1, 0.8)

  ord <- order(time)
  time_s <- time[ord]
  X_s <- cbind(x1, x2)[ord, ]
  st_s <- status[ord]

  b <- cbind(rep(0, 2))
  for(i in 1:50) {
    eta <- c(X_s %**% b)
    sc <- score_cox(X_s, eta, st_s, time_s)
    info <- info_cox(X_s, eta, st_s, time_s)
    step <- solve(info, sc)
    b <- b + c(step)
  }

  eta_final <- c(X_s %**% b)
  sc_final <- score_cox(X_s, eta_final, st_s, time_s)
  expect_true(all(abs(sc_final) < 1e-4))

  ## Test 3: Compare with survival::coxph
  skip_if_not_installed("survival")
  df_test <- data.frame(time = time, status = status, x1 = x1, x2 = x2)
  coxfit <- survival::coxph(
    survival::Surv(time, status) ~ x1 + x2,
    data = df_test, ties = "breslow"
  )
  b_coxph <- coef(coxfit)
  names(b_coxph) <- NULL
  expect_equal(c(b), b_coxph, tolerance = 1e-2)

  ## Test 4: Information matrix is positive definite at the MLE
  eta_mle <- c(X_s %**% b)
  I <- info_cox(X_s, eta_mle, st_s, time_s)
  eigvals <- eigen(I, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigvals > 0))

  ## Test 5: Partial log-likelihood increases along Newton direction
  b0 <- cbind(rep(0, 2))
  eta0 <- c(X_s %**% b0)
  ll0 <- loglik_cox(eta0, st_s, time_s)

  sc0 <- score_cox(X_s, eta0, st_s, time_s)
  I0  <- info_cox(X_s, eta0, st_s, time_s)
  step0 <- solve(I0, sc0)
  b1 <- b0 + c(step0)
  eta1 <- c(X_s %**% b1)
  ll1 <- loglik_cox(eta1, st_s, time_s)
  expect_gt(ll1, ll0)

  ## Test 6: unconstrained_fit_cox recovers MLE (no penalty)
  p <- 2
  Lambda_zero <- matrix(0, p, p)
  LambdaHalf_zero <- matrix(0, p, p)

  b_fit <- unconstrained_fit_cox(
    X = cbind(x1, x2), y = time,
    LambdaHalf = LambdaHalf_zero, Lambda = Lambda_zero,
    keep_weighted_Lambda = FALSE, family = cox_family(),
    tol = 1e-10, K = 0,
    parallel = FALSE, cl = NULL,
    chunk_size = NULL, num_chunks = NULL, rem_chunks = NULL,
    order_indices = 1:n, weights = rep(1, n), status = status
  )
  expect_equal(c(b_fit), c(b), tolerance = 1e-2)

  ## Test 7: cox_qp_score_function near zero at MLE
  X_full <- cbind(x1, x2)
  mu_test <- exp(c(X_full %**% b))
  order_list_test <- list(1:n)

  sc_qp <- cox_qp_score_function(
    X = X_full, y = time, mu = mu_test,
    order_list = order_list_test, dispersion = 1,
    VhalfInv = NULL, observation_weights = rep(1, n),
    status = status
  )
  expect_true(all(abs(sc_qp) < 1e-4))

  ## Test 8: Weighted partial log-likelihood
  set.seed(99)
  n_wt <- 20
  time_wt <- rexp(n_wt, 1)
  eta_wt <- rnorm(n_wt, 0, 0.3)
  st_wt <- rbinom(n_wt, 1, 0.6)
  wt_wt <- runif(n_wt, 0.5, 2.0)

  ord_wt <- order(time_wt)
  time_wt_s <- time_wt[ord_wt]
  eta_wt_s <- eta_wt[ord_wt]
  st_wt_s <- st_wt[ord_wt]
  wt_wt_s <- wt_wt[ord_wt]

  haz_wt <- wt_wt_s * exp(eta_wt_s)
  event_times_wt <- sort(unique(time_wt_s[st_wt_s == 1]))
  ll_manual <- 0
  for(tt in event_times_wt) {
    Dg <- which(time_wt_s == tt & st_wt_s == 1)
    Rg <- which(time_wt_s >= tt)
    dgw <- sum(wt_wt_s[Dg])
    ll_manual <- ll_manual +
      sum(wt_wt_s[Dg] * eta_wt_s[Dg]) -
      dgw * log(sum(haz_wt[Rg]))
  }
  ll_fn <- loglik_cox(eta_wt_s, st_wt_s, time_wt_s, wt_wt_s)
  expect_equal(ll_fn, ll_manual, tolerance = 1e-4)

  ll_unit <- loglik_cox(eta_wt_s, st_wt_s, time_wt_s, rep(1, n_wt))
  ll_default <- loglik_cox(eta_wt_s, st_wt_s, time_wt_s)
  expect_equal(ll_unit, ll_default, tolerance = 1e-4)

  ## Test 9: Dispersion function returns 1
  d_val <- cox_dispersion_function(mu_test, time, 1:n, cox_family(),
                                   rep(1, n), NULL)
  expect_identical(d_val, 1)

  ## Test 10: GLM weight function returns positive values
  wt_vals <- cox_glm_weight_function(
    mu = mu_test, y = time, order_indices = 1:n,
    family = cox_family(), dispersion = 1,
    observation_weights = rep(1, n), status = status
  )
  expect_true(all(wt_vals > 0))
  expect_length(wt_vals, n)
})
