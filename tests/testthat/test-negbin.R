test_that("Negative Binomial helpers", {
  ## Test 1: Log-likelihood at known values matches manual computation
  set.seed(1234)
  N <- 50
  y <- rpois(N, 5)
  mu <- rep(5, N)
  theta <- 3
  ll_manual <- sum(
    lgamma(y + theta) - lgamma(theta) - lgamma(y + 1) +
      theta * log(theta) - theta * log(mu + theta) +
      y * log(mu) - y * log(mu + theta)
  )
  ll_fn <- loglik_negbin(y, mu, theta)
  expect_equal(ll_fn, ll_manual, tolerance = 1e-4)

  ## Test 2: Score at MLE should be approximately zero
  set.seed(42)
  N <- 300
  t1 <- rnorm(N)
  t2 <- rbinom(N, 1, 0.5)
  true_theta <- 4
  mu_true <- exp(1 + 0.3 * t1 + 0.2 * t2)
  y <- rnbinom(N, size = true_theta, mu = mu_true)
  X_s <- cbind(1, t1, t2)
  b <- cbind(c(log(mean(y)), 0, 0))
  ## Profile: alternate theta and beta
  for(outer in 1:10) {
    mu_cur <- c(exp(X_s %**% b))
    th <- negbin_theta(y, mu_cur)
    for(i in 1:20) {
      mu_cur <- c(exp(X_s %**% b))
      sc <- score_negbin(X_s, y, mu_cur, th)
      I <- info_negbin(X_s, y, mu_cur, th)
      step <- solve(I, sc)
      b <- b + c(step)
    }
  }
  mu_final <- c(exp(X_s %**% b))
  th_final <- negbin_theta(y, mu_final)
  sc_final <- score_negbin(X_s, y, mu_final, th_final)
  expect_true(all(abs(sc_final) < 1e-4))

  ## Test 3: Compare with MASS::glm.nb
  skip_if_not_installed("MASS")
  df_test <- data.frame(y = y, t1 = t1, t2 = t2)
  nbfit <- MASS::glm.nb(y ~ t1 + t2, data = df_test)
  b_mass <- coef(nbfit)
  names(b_mass) <- NULL
  expect_equal(c(b), b_mass, tolerance = 5e-2)
  ## Theta comparison
  expect_equal(th_final, nbfit$theta, tolerance = 1)

  ## Test 4: Information matrix is positive definite at the MLE
  I_mle <- info_negbin(X_s, y, mu_final, th_final)
  eigvals <- eigen(I_mle, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigvals > 0))

  ## Test 5: Log-likelihood increases along Newton direction
  b0 <- cbind(c(log(mean(y)), 0, 0))
  mu0 <- c(exp(X_s %**% b0))
  th0 <- negbin_theta(y, mu0)
  ll0 <- loglik_negbin(y, mu0, th0)
  sc0 <- score_negbin(X_s, y, mu0, th0)
  I0 <- info_negbin(X_s, y, mu0, th0)
  step0 <- solve(I0, sc0)
  b1 <- b0 + c(step0)
  mu1 <- c(exp(X_s %**% b1))
  th1 <- negbin_theta(y, mu1)
  ll1 <- loglik_negbin(y, mu1, th1)
  expect_gt(ll1, ll0)

  ## Test 6: unconstrained_fit_negbin recovers MLE (no penalty)
  p <- 3
  Lambda_zero <- matrix(0, p, p)
  LambdaHalf_zero <- matrix(0, p, p)
  b_fit <- unconstrained_fit_negbin(
    X = X_s, y = y,
    LambdaHalf = LambdaHalf_zero, Lambda = Lambda_zero,
    keep_weighted_Lambda = FALSE, family = negbin_family(),
    tol = 1e-10, K = 0,
    parallel = FALSE, cl = NULL,
    chunk_size = NULL, num_chunks = NULL, rem_chunks = NULL,
    order_indices = 1:N, weights = rep(1, N)
  )
  expect_equal(c(b_fit), c(b), tolerance = 5e-2)

  ## Test 7: negbin_qp_score_function near zero at MLE
  mu_test <- exp(c(X_s %**% b))
  order_list_test <- list(1:N)
  sc_qp <- negbin_qp_score_function(
    X = X_s, y = y, mu = mu_test,
    order_list = order_list_test, dispersion = th_final,
    VhalfInv = NULL, observation_weights = rep(1, N)
  )
  expect_true(all(abs(sc_qp) < 1e-4))

  ## Test 8: Weighted log-likelihood
  set.seed(99)
  N_wt <- 30
  y_wt <- rpois(N_wt, 3)
  mu_wt <- rep(3, N_wt)
  th_wt <- 5
  wt_wt <- runif(N_wt, 0.5, 2.0)
  ll_manual <- sum(wt_wt * (
    lgamma(y_wt + th_wt) - lgamma(th_wt) - lgamma(y_wt + 1) +
      th_wt * log(th_wt) - th_wt * log(mu_wt + th_wt) +
      y_wt * log(mu_wt) - y_wt * log(mu_wt + th_wt)
  ))
  ll_fn <- loglik_negbin(y_wt, mu_wt, th_wt, wt_wt)
  expect_equal(ll_fn, ll_manual, tolerance = 1e-4)
  ll_unit <- loglik_negbin(y_wt, mu_wt, th_wt, rep(1, N_wt))
  ll_default <- loglik_negbin(y_wt, mu_wt, th_wt)
  expect_equal(ll_unit, ll_default, tolerance = 1e-4)

  ## Test 9: Dispersion function returns positive theta
  th_est <- negbin_dispersion_function(
    mu = mu_test, y = y, order_indices = 1:N,
    family = negbin_family(),
    observation_weights = rep(1, N),
    VhalfInv = NULL
  )
  expect_true(th_est > 0)
  expect_true(is.finite(th_est))

  ## Test 10: GLM weight function returns positive values
  wt_vals <- negbin_glm_weight_function(
    mu = mu_test, y = y, order_indices = 1:N,
    family = negbin_family(), dispersion = th_final,
    observation_weights = rep(1, N)
  )
  expect_true(all(wt_vals > 0))
  expect_length(wt_vals, N)

  ## Test 11: negbin_theta recovers known theta from large sample
  set.seed(2025)
  N_big <- 5000
  mu_big <- rep(10, N_big)
  theta_true <- 3.0
  y_big <- rnbinom(N_big, size = theta_true, mu = 10)
  th_hat <- negbin_theta(y_big, mu_big)
  expect_equal(th_hat, theta_true, tolerance = 0.5)

  ## Test 12: Schur correction is negative semi-definite
  X_test <- cbind(1, t1, t2)
  B_test <- list(b)
  order_list_sc <- list(1:N)
  obs_wt_sc <- list(rep(1, N))
  schur <- negbin_schur_correction(
    X = list(X_test), y = list(y),
    B = B_test, dispersion = th_final,
    order_list = order_list_sc, K = 0,
    family = negbin_family(),
    observation_weights = obs_wt_sc
  )
  ## Should be a matrix (not 0) for nondegenerate data
  expect_true(is.matrix(schur[[1]]))
  ## Eigenvalues of Schur correction should be <= 0 (negative semi-definite)
  eig_schur <- eigen(schur[[1]], symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig_schur <= sqrt(.Machine$double.eps)))

  ## Test 13: Score with weights is consistent
  set.seed(77)
  wt_sc <- runif(N, 0.5, 2.0)
  mu_sc <- c(exp(X_s %**% b))
  sc_weighted <- score_negbin(X_s, y, mu_sc, th_final, wt_sc)
  ## Manual: X^T (wt * (y - mu) * theta / (theta + mu))
  r_manual <- wt_sc * (y - mu_sc) * th_final / (th_final + mu_sc)
  sc_manual <- crossprod(X_s, cbind(r_manual))
  expect_equal(c(sc_weighted), c(sc_manual), tolerance = 1e-4)

  ## Test 14: Info with weights matches X^T W X manually
  W_manual <- wt_sc * mu_sc * th_final * (th_final + y) / (th_final + mu_sc)^2
  I_manual <- crossprod(X_s, W_manual * X_s)
  I_fn <- info_negbin(X_s, y, mu_sc, th_final, wt_sc)
  expect_equal(I_fn, I_manual, tolerance = 1e-4)

  ## Test 15: Dispersion function with VhalfInv gives valid result
  N_gee <- 50
  y_gee <- rnbinom(N_gee, size = 3, mu = 5)
  mu_gee <- rep(5, N_gee)
  rho <- 0.2
  ## Simple diagonal perturbation for testing
  VhalfInv_test <- diag(1 / sqrt(1 - rho), N_gee, N_gee)
  th_gee <- negbin_dispersion_function(
    mu = mu_gee, y = y_gee, order_indices = 1:N_gee,
    family = negbin_family(),
    observation_weights = rep(1, N_gee),
    VhalfInv = VhalfInv_test
  )
  expect_true(th_gee > 0)
  expect_true(is.finite(th_gee))

  ## Test 16: Poisson limit (theta -> Inf gives Poisson log-likelihood)
  set.seed(88)
  N_pois <- 100
  mu_pois <- rep(3, N_pois)
  y_pois <- rpois(N_pois, 3)
  ll_nb_big <- loglik_negbin(y_pois, mu_pois, theta = 1e8)
  ll_pois <- sum(dpois(y_pois, mu_pois, log = TRUE))
  expect_equal(ll_nb_big, ll_pois, tolerance = 1e-2)

  ## Test 17: keep_weighted_Lambda returns Poisson hot-start
  b_hot <- unconstrained_fit_negbin(
    X = X_s, y = y,
    LambdaHalf = diag(0.01, p), Lambda = diag(0.01^2, p),
    keep_weighted_Lambda = TRUE, family = negbin_family(),
    tol = 1e-8, K = 0,
    parallel = FALSE, cl = NULL,
    chunk_size = NULL, num_chunks = NULL, rem_chunks = NULL,
    order_indices = 1:N, weights = rep(1, N)
  )
  expect_length(b_hot, p)
  expect_true(all(is.finite(b_hot)))

  ## Test 18: Family object components exist and are correct
  fam <- negbin_family()
  expect_equal(fam$family, "negbin")
  expect_equal(fam$link, "log")
  expect_equal(fam$linkfun(exp(1)), 1)
  expect_equal(fam$linkinv(0), 1)
  expect_true(is.function(fam$loglik))
  expect_true(is.function(fam$aic))
  expect_true(is.function(fam$dev.resids))
  expect_true(is.function(fam$custom_dev.resids))

  ## Test 19: Schur correction returns 0 for empty partition
  schur_empty <- negbin_schur_correction(
    X = list(matrix(0, nrow = 0, ncol = 3), X_test),
    y = list(numeric(0), y),
    B = list(cbind(rep(0, 3)), b),
    dispersion = th_final,
    order_list = list(integer(0), 1:N),
    K = 1,
    family = negbin_family(),
    observation_weights = list(numeric(0), rep(1, N))
  )
  expect_identical(schur_empty[[1]], 0)
  expect_true(is.matrix(schur_empty[[2]]))

  ## Test 20: Numerical gradient check for score_negbin
  set.seed(123)
  X_ng <- cbind(1, rnorm(80))
  y_ng <- rnbinom(80, size = 2, mu = 3)
  b_ng <- c(log(3), 0.1)
  th_ng <- 2
  eps_ng <- 1e-6
  sc_analytic <- score_negbin(X_ng, y_ng, exp(X_ng %**% cbind(b_ng)), th_ng)
  sc_numerical <- sapply(1:2, function(j) {
    bp <- bm <- b_ng
    bp[j] <- bp[j] + eps_ng
    bm[j] <- bm[j] - eps_ng
    (loglik_negbin(y_ng, exp(X_ng %**% cbind(bp)), th_ng) -
        loglik_negbin(y_ng, exp(X_ng %**% cbind(bm)), th_ng)) / (2 * eps_ng)
  })
  expect_equal(c(sc_analytic), sc_numerical, tolerance = 1e-4)
})
