test_that("1D complex function examples run without error", {
  expect_no_error({

    ## 1-D Setup
    set.seed(1234)
    x <- seq(-9, 9, length.out = 1000)
    slinky <- function(x) (50 * cos(x * 2) +-2 * x ^ 2 + (0.25 * x) ^ 4 + 80)
    coil <- function(x) (100 * cos(x * 2) +-1.5 * x ^ 2 + (0.1 * x) ^ 4 + (0.05 * x ^ 3) +
                           (-0.01 * x ^ 5) + (0.00002*x^6) -(0.000001*x^7) + 100)
    exponential_log <- function(x) {
      unlist(c(sapply(x, function(xx) {
        if (xx <= 1) {
          100 * (exp(xx) - exp(1))
        } else {
          100 * (log(xx))
        }
      })))
    }
    scaled_abs_gamma <- function(x) 2*sqrt(gamma(abs(x)))
    fxn <- function(x)(slinky(x) + coil(x) + exponential_log(x) + scaled_abs_gamma(x))
    dat <- cbind(x, fxn(x) + rnorm(length(x), 0, 50))
    colnames(dat) <- c('x', 'y')


    ## 3-ways to fit model, 3rd is the basic functionality without
    # plotting, posterior generation, wald inference, or extremum finding
    fit1 <- lgspline(predictors = cbind(dat[,'x']),
                     use_custom_bfgs = FALSE,
                     response = dat[,'y'])
    fit2 <- lgspline(y ~ spl(x),
                     data.frame(dat))
    fit3 <- lgspline.fit(cbind(dat[,'x']),
                         dat[,'y'])

    ## Plotting
    plot(fit2)
    points(dat)
    fit1$plot(show_formulas = TRUE,
              digits = 1,
              new_predictors = cbind(seq(-10, 10, length.out = 10000)),
              ylim = c(-1000, 7500),
              xlim = c(-11, 11),
              legend_pos = 'topleft',
              custom_response_lab = 'Batman',
              custom_predictor_lab = 'Robin',
              custom_formula_lab = 'TheJoker',
              custom_title = 'Playing with some non-default settings',
              text_size_formula = 0.5)

    ## No first derivative constraints, 1 knot, constrain intercept = 0,
    # include weights
    fit4 <- lgspline(y ~ 0 + spl(x),
                     standardize_response = FALSE,
                     data.frame(dat),
                     opt = FALSE,
                     weights = abs(rnorm(nrow(dat)))*0.5 + 0.5,
                     K = 1,
                     include_constrain_first_deriv = FALSE)
    fit4$plot(custom_title =
                'Exclude first derivative constraints and intercept',
            show_formulas = TRUE)

    ## Leave-one-out cross-validated predictions
    plot(leave_one_out(fit1), fit1$y)
    abline(0, 1)
  })
})


test_that("trees data examples run without error", {
  expect_no_error({

    ## 2-D Setup
    set.seed(1234)
    data("trees")

    ## Test expansions_only argument prevents
    exp_only <- lgspline(
      Volume ~ spl(Girth) + Height*Girth,
      trees,
      log_initial_flat = 1,
      log_initial_wiggle = 1e-8,
      K = 1,
      expansions_only = TRUE,
      family = Gamma(link = 'log')
    )
    expect_null(exp_only$B)

    ## Fit model assuming gamma-distributed outcome
    model_fit <- lgspline(
      Volume ~ spl(Girth) + Height*Girth,
      trees,
      log_initial_flat = 1,
      log_initial_wiggle = 1e-8,
      K = 1,
      opt = FALSE,
      family = Gamma(link = 'log'),
      custom_knots = exp_only$knots,
      make_partition_list = exp_only$make_partition_list
    )

    ## Some basic functionality
    predict(model_fit)
    print(summary(model_fit))
    generate_posterior(model_fit)
    find_extremum(model_fit)
    model_fit$plot()

    ## Use quadratic-programming to enforce bounds on absolute value of
    # all coefficients
    # Custom constraint matrix
    l1_constraint_matrix <- function(p, K) {
      P <- p * (K + 1)
      first_diag <- diag(P)
      second_diag <- -diag(P)
      l1_Amat <- cbind(first_diag, second_diag)
      return(l1_Amat)
    }
    # Custom constraint bounds
    l1_bound_vector <- function(qp_Amat, scales, l1_bound) {
      l1_bvec <- rep(-l1_bound, ncol(qp_Amat)) * c(1, scales)
      return(l1_bvec)
    }

    ## Fit the same model but with some additional custom options
    # 4 partitions
    l1_bound <- 0.1
    model_fit <- lgspline(
      predictors = with(trees, cbind(Girth, Height)),
      response = trees$Volume,
      family = Gamma(link = 'log'),
      just_linear_with_interactions = 1,
      neighbor_tolerance = 2,
      wiggle_penalty = 1,
      flat_ridge_penalty = 1,
      unique_penalty_per_partition = FALSE,
      unique_penalty_per_predictor = FALSE,
      opt = FALSE,
      K = 4,
      qp_Amat_fxn = function(N, p, K, X_block, colnm, scales, fxn, ...) {
        mat <- l1_constraint_matrix(p, K)
        mat
      },
      qp_bvec_fxn = function(qp_Amat, N, p, K, X_block, colnm, scales, fxn, ...) {
        vec <- l1_bound_vector(qp_Amat, scales, l1_bound)
        vec
      },
      qp_meq_fxn = function(qp_Amat, N, p, K, X_block, colnm, scales, fxn, ... ) 0
    )

    ## Plotting with some custom options
    model_fit$plot(custom_predictor_lab1 = 'Girth',
                   custom_predictor_lab2 = 'Height',
                   custom_response_lab = 'Volume',
                   custom_title = 'Girth and Height Predicting Volume of Trees',
                   show_formulas = TRUE,
                   new_predictors = expand.grid(
                     seq(min(trees$Girth), max(trees$Girth), length.out = 250),
                     seq(min(trees$Height), max(trees$Height), length.out = 250)
                   ))

    ## Check the absolute magnitude of all coefficients do not surpass 0.1
    if(max(abs(unlist(model_fit$B))) > 0.1+1e-8) stop('Inequality not met')

    ## Find the predictor yielding fitted values most close to the median,
    # using a custom objective function and gradient
    find_extremum(
      model_fit,
      minimize = TRUE,
      custom_objective_function = function(mu, sigma, ybest){
        0.5*(mu - median(trees$Volume))^2
      },
      custom_objective_derivative = function(mu, sigma, ybest, d_mu){
        (mu - median(trees$Volume)) * d_mu
      }
    )

  })
})

test_that("simple linear regression approximately matches base R inference", {
  set.seed(1234)

  ## Running simple linear regression, outputs should approximately match
  x <- rnorm(100)
  y <- x + rnorm(100)
  dat <- data.frame(
    x = x,
    y = y
  )
  fit_lgspline <- lgspline(formula = y ~ 0 + x,
                           wiggle_penalty = 0,
                           flat_ridge_penalty = 0,
                           opt = FALSE,
                           unique_penalty_per_partition = FALSE,
                           unique_penalty_per_predictor = FALSE,
                           data = data.frame(dat))
  coefs <- c(unlist(coef(fit_lgspline)))
  names(coefs) <- NULL
  lms <- summary(lm(y ~ 0 + x, data = data.frame(dat)))

  ## Estimates should strongly match
  expect_equal(round(c(0, lms$coefficients[1]), 10),
               round(coefs, 10),
               tolerance = 1e-4)

  ## So should dispersion estimates, but to a lesser extent
  expect_equal(round(lms$sigma, 10),
               round(sqrt(fit_lgspline$sigmasq_tilde), 10),
               tolerance = 1e-2)

  ## Ok for slightly-different standard errors
  expect_equal(round(lms$coefficients[2], 10),
               round(fit_lgspline$wald_univariate()$se[2], 10),
               tolerance = 1e-2)

})

test_that("Weibull AFT models can be fit under custom constraints", {
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
  model_fit <- lgspline(y ~ spl(x1) + x2,
                        df,
                        unconstrained_fit_fxn = unconstrained_fit_weibull,
                        family = weibull_family(),
                        need_dispersion_for_estimation = TRUE,
                        qp_score_function = weibull_qp_score_function(),
                        dispersion_function = weibull_dispersion_function,
                        glm_weight_function = weibull_glm_weight_function,
                        shur_correction_function = weibull_shur_correction,
                        K = 1,
                        opt = FALSE,
                        return_varcovmat = TRUE,
                        observation_weights = abs(rnorm(1000)),
                        constraint_vectors = cbind(rep(1, 12)),
                        null_constraint = matrix(12),
                        status = status,
                        verbose = TRUE)

  ## Check sum-to-P constraint
  expect_equal(round(sum(unlist(model_fit$B)), 10), 12)
})
