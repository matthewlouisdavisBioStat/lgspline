set.seed(1234)
library(lgspline)

################################################################################

## ## 95% Coverage Operating Characteristics ## ##
#  This script runs simulation checks for coefficient and correlation-parameter
#  coverage across the main lgspline likelihood paths. It is intentionally not a
#  CRAN test: at reviewer-scale settings it can take a long time.
#
#  The common regression target is
#
#      y ~ spl(t1, t2) + t1*t3
#
#  so the spline surface is a joined two-predictor smooth and the coefficient
#  target is the explicit t1-by-t3 interaction, repeated over all partitions.
#
#  Environment overrides:
#    LGSPLINE_COVERAGE_N          sample size; default 2000
#    LGSPLINE_COVERAGE_REPS       Monte Carlo replicates; default 100
#    LGSPLINE_COVERAGE_FAST       TRUE uses K = 1 and opt = FALSE for smoke runs
#    LGSPLINE_COVERAGE_FAST_K     K used in fast smoke mode; default 1
#    LGSPLINE_COVERAGE_SCENARIOS  comma-separated scenario names
#    LGSPLINE_COVERAGE_CSV        output CSV path

## Simulation controls
N <- as.integer(Sys.getenv("LGSPLINE_COVERAGE_N", "2000"))
n_sims <- as.integer(Sys.getenv("LGSPLINE_COVERAGE_REPS", "100"))
fast_smoke <- tolower(Sys.getenv("LGSPLINE_COVERAGE_FAST", "FALSE")) %in%
  c("true", "t", "1", "yes", "y")
fast_K <- as.integer(Sys.getenv("LGSPLINE_COVERAGE_FAST_K", "1"))
scenario_filter <- Sys.getenv("LGSPLINE_COVERAGE_SCENARIOS", "")
output_csv <- Sys.getenv(
  "LGSPLINE_COVERAGE_CSV",
  file.path(getwd(), "coverage_operating_characteristics.csv")
)

## True parameter values
beta_int <- 0.35
beta_t3 <- -0.25
sigma_gaussian <- 0.75
gamma_shape <- 6
negbin_size <- 4
weibull_scale_true <- 0.75

## Correlation parameters are stored by lgspline on working scales.
rho_ar1 <- 0.45
rho_exchangeable <- 0.25
omega_spatial <- 1.10
matern_ell <- 0.70
matern_nu <- 1.25
true_working <- list(
  ar1 = log(-log(rho_ar1)),
  exchangeable = log(-log(rho_exchangeable)),
  spatial_exponential = log(omega_spatial),
  matern = log(c(matern_ell, matern_nu))
)

## Fit controls: leave defaults alone unless doing a smoke run.
base_fit_args <- list(
  return_varcovmat = TRUE,
  include_warnings = FALSE
)
if(fast_smoke) {
  base_fit_args$K <- fast_K
  base_fit_args$opt <- FALSE
}

################################################################################

## ## Data-Generating Helpers ## ##

make_covariates <- function(N, block_size = NULL) {
  dat <- data.frame(
    t1 = runif(N, -1, 1),
    t2 = runif(N, -1, 1),
    t3 = runif(N, -1, 1)
  )

  if(!is.null(block_size)) {
    if(N %% block_size != 0) {
      stop("N must be divisible by block_size for correlated scenarios.")
    }
    dat$id <- rep(seq_len(N / block_size), each = block_size)
    dat$spacetime <- rep(seq(0, 1, length.out = block_size), N / block_size)
  }
  dat
}

smooth_signal <- function(t1, t2) {
  0.35 * sin(2.1 * t1 + 0.3 * cos(3 * t2)) +
    0.25 * cos(1.7 * t2 - 0.4 * t1^2) -
    0.30 * exp(-((t1 - 0.35)^2 + (t2 + 0.25)^2) / 0.18) +
    0.12 * atan(3 * t1 * t2)
}

eta_common <- function(dat, intercept = 0, scale = 1) {
  intercept +
    scale * smooth_signal(dat$t1, dat$t2) +
    beta_t3 * dat$t3 +
    beta_int * dat$t1 * dat$t3
}

block_cor_ar1 <- function(d, rho = rho_ar1) {
  ranked <- abs(outer(seq_along(d), seq_along(d), "-"))
  rho^ranked
}

block_cor_exchangeable <- function(d, rho = rho_exchangeable) {
  m <- length(d)
  (1 - rho) * diag(m) + rho * matrix(1, m, m)
}

block_cor_spatial_exponential <- function(d, omega = omega_spatial) {
  exp(-omega * abs(outer(d, d, "-")))
}

block_cor_matern <- function(d, ell = matern_ell, nu = matern_nu) {
  diffs <- abs(outer(d, d, "-"))
  scaled <- sqrt(2 * nu) * diffs / ell
  V <- matrix(1, length(d), length(d))
  nz <- scaled > 0
  V[nz] <- (2^(1 - nu) / gamma(nu)) *
    scaled[nz]^nu * besselK(scaled[nz], nu)
  V
}

make_exchangeable_vhalf_functions <- function(id) {
  list(
    VhalfInv_fxn = function(par) {
      rho <- exp(-exp(par))
      out <- matrix(0, length(id), length(id))
      for(clust in unique(id)) {
        inds <- which(id == clust)
        m <- length(inds)
        V <- (1 - rho) * diag(m) + rho * matrix(1, m, m)
        out[inds, inds] <- matinvsqrt(V)
      }
      out
    },
    Vhalf_fxn = function(par) {
      rho <- exp(-exp(par))
      out <- matrix(0, length(id), length(id))
      for(clust in unique(id)) {
        inds <- which(id == clust)
        m <- length(inds)
        V <- (1 - rho) * diag(m) + rho * matrix(1, m, m)
        out[inds, inds] <- matsqrt(V)
      }
      out
    },
    VhalfInv_logdet = function(par) {
      rho <- exp(-exp(par))
      log_det <- 0
      for(clust in unique(id)) {
        inds <- which(id == clust)
        m <- length(inds)
        V <- (1 - rho) * diag(m) + rho * matrix(1, m, m)
        log_det <- log_det -
          0.5 * determinant(V, logarithm = TRUE)$modulus[1]
      }
      log_det
    }
  )
}

make_spatial_exponential_vhalf_functions <- function(id, spacetime) {
  list(
    VhalfInv_fxn = function(par) {
      omega <- exp(par)
      out <- matrix(0, length(id), length(id))
      for(clust in unique(id)) {
        inds <- which(id == clust)
        d <- abs(outer(spacetime[inds], spacetime[inds], "-"))
        out[inds, inds] <- matinvsqrt(exp(-omega * d))
      }
      out
    },
    Vhalf_fxn = function(par) {
      omega <- exp(par)
      out <- matrix(0, length(id), length(id))
      for(clust in unique(id)) {
        inds <- which(id == clust)
        d <- abs(outer(spacetime[inds], spacetime[inds], "-"))
        out[inds, inds] <- matsqrt(exp(-omega * d))
      }
      out
    },
    VhalfInv_logdet = function(par) {
      omega <- exp(par)
      log_det <- 0
      for(clust in unique(id)) {
        inds <- which(id == clust)
        d <- abs(outer(spacetime[inds], spacetime[inds], "-"))
        V <- exp(-omega * d)
        log_det <- log_det -
          0.5 * determinant(V, logarithm = TRUE)$modulus[1]
      }
      log_det
    }
  )
}

draw_block_normal <- function(id, spacetime, block_cor_fxn) {
  z <- numeric(length(id))
  for(clust in unique(id)) {
    inds <- which(id == clust)
    V <- block_cor_fxn(spacetime[inds])
    R <- chol(V)
    z[inds] <- c(t(R) %*% rnorm(length(inds)))
  }
  z
}

draw_correlated_count <- function(mu, id, spacetime, block_cor_fxn) {
  ## Count-valued quasi-Poisson DGP with approximately correct working
  #  Pearson-residual correlation:
  #      y_i = round(mu_i + sqrt(mu_i) * z_i), z ~ N(0, V).
  #  The pmax() truncation is negligible for the moderately large means used
  #  below, but keeps the response in the count support.
  z <- draw_block_normal(id, spacetime, block_cor_fxn)
  as.integer(pmax(0, round(mu + sqrt(pmax(mu, .Machine$double.eps)) * z)))
}

simulate_weibull_aft <- function(mu, scale) {
  ## log(T) = log(mu) + scale * log(E), E ~ Exp(1)
  mu * rexp(length(mu))^scale
}

apply_censoring <- function(time, rate) {
  cens <- rexp(length(time), rate = rate)
  list(
    y = pmin(time, cens),
    status = as.integer(time <= cens)
  )
}

################################################################################

## ## Fit and Extraction Helpers ## ##

call_fit <- function(fit_fun, formula, data, extra_args = list()) {
  args <- c(list(formula = formula, data = data), base_fit_args, extra_args)
  do.call(fit_fun, args)
}

extract_interaction_coverage <- function(fit, scenario, sim_id) {
  ci <- confint(fit)
  wald <- wald_univariate(fit)$coefficients
  rows <- grep("(^|\\.)t1xt3$|(^|\\.)t3xt1$", rownames(ci), value = TRUE)
  if(length(rows) == 0) {
    rows <- grep("t1xt3|t3xt1", rownames(ci), value = TRUE)
  }

  if(length(rows) == 0) {
    return(data.frame(
      scenario = scenario,
      sim = sim_id,
      estimand = "interaction",
      parameter = "t1xt3",
      true = beta_int,
      estimate = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      covered = NA,
      status = "missing_interaction_ci"
    ))
  }

  lower <- ci[rows, 1]
  upper <- ci[rows, 2]
  covered <- lower <= beta_int & beta_int <= upper
  status <- ifelse(is.finite(lower) & is.finite(upper),
                   "ok", "nonfinite_interaction_ci")

  data.frame(
    scenario = scenario,
    sim = sim_id,
    estimand = "interaction",
    parameter = rows,
    true = beta_int,
    estimate = wald[rows, "Estimate"],
    lower = lower,
    upper = upper,
    covered = covered,
    status = status,
    row.names = NULL
  )
}

extract_correlation_coverage <- function(fit, scenario, sim_id, true_par) {
  ci <- confint(fit)
  rows <- grep("^Correlation parameter", rownames(ci), value = TRUE)
  if(length(rows) == 0) {
    return(data.frame(
      scenario = scenario,
      sim = sim_id,
      estimand = "correlation",
      parameter = paste0("Correlation parameter ", seq_along(true_par)),
      true = true_par,
      estimate = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      covered = NA,
      status = "missing_correlation_ci"
    ))
  }

  true_par <- rep(true_par, length.out = length(rows))
  lower <- ci[rows, 1]
  upper <- ci[rows, 2]
  covered <- lower <= true_par & true_par <= upper
  status <- ifelse(is.finite(lower) & is.finite(upper),
                   "ok", "nonfinite_correlation_ci")

  data.frame(
    scenario = scenario,
    sim = sim_id,
    estimand = "correlation",
    parameter = rows,
    true = true_par,
    estimate = fit$VhalfInv_params_estimates[seq_along(rows)],
    lower = lower,
    upper = upper,
    covered = covered,
    status = status,
    row.names = NULL
  )
}

summarize_results <- function(results) {
  ok <- results[results$status == "ok" & !is.na(results$covered), , drop = FALSE]
  if(nrow(ok) == 0) return(data.frame())

  split_res <- split(ok, list(ok$scenario, ok$estimand, ok$parameter),
                     drop = TRUE)
  out <- do.call(rbind, lapply(split_res, function(x) {
    p <- mean(x$covered)
    data.frame(
      scenario = x$scenario[1],
      estimand = x$estimand[1],
      parameter = x$parameter[1],
      n = nrow(x),
      coverage = p,
      mc_se = sqrt(p * (1 - p) / nrow(x)),
      mean_estimate = mean(x$estimate),
      mean_width = mean(x$upper - x$lower),
      row.names = NULL
    )
  }))
  out[order(out$scenario, out$estimand, out$parameter), ]
}

write_results <- function(results, output_csv) {
  if(is.null(results) || nrow(results) == 0) return(invisible(NULL))
  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  write.csv(results, output_csv, row.names = FALSE)
  summary_csv <- sub("\\.csv$", "_summary.csv", output_csv)
  write.csv(summarize_results(results), summary_csv, row.names = FALSE)
  invisible(summary_csv)
}

record_failure <- function(scenario, sim_id, err) {
  data.frame(
    scenario = scenario,
    sim = sim_id,
    estimand = "fit",
    parameter = "fit",
    true = NA_real_,
    estimate = NA_real_,
    lower = NA_real_,
    upper = NA_real_,
    covered = NA,
    status = paste("error:", conditionMessage(err)),
    row.names = NULL
  )
}

################################################################################

## ## Scenario Definitions ## ##

fit_formula <- y ~ spl(t1, t2) + t1*t3
time_formula <- time ~ spl(t1, t2) + t1*t3

scenario_gaussian <- function(sim_id) {
  dat <- make_covariates(N)
  eta <- eta_common(dat, intercept = 0.2, scale = 1)
  dat$y <- eta + rnorm(N, 0, sigma_gaussian)
  fit <- call_fit(lgspline, fit_formula, dat)
  extract_interaction_coverage(fit, "gaussian", sim_id)
}

scenario_logistic <- function(sim_id) {
  dat <- make_covariates(N)
  eta <- eta_common(dat, intercept = -0.10, scale = 0.75)
  dat$y <- rbinom(N, 1, plogis(eta))
  fit <- call_fit(lgspline, fit_formula, dat,
                  list(family = binomial()))
  extract_interaction_coverage(fit, "logistic", sim_id)
}

scenario_poisson <- function(sim_id) {
  dat <- make_covariates(N)
  eta <- eta_common(dat, intercept = 0.55, scale = 0.65)
  dat$y <- rpois(N, exp(eta))
  fit <- call_fit(lgspline, fit_formula, dat,
                  list(family = poisson()))
  extract_interaction_coverage(fit, "poisson", sim_id)
}

scenario_gamma_log <- function(sim_id) {
  dat <- make_covariates(N)
  eta <- eta_common(dat, intercept = 0.45, scale = 0.55)
  mu <- exp(eta)
  dat$y <- rgamma(N, shape = gamma_shape, scale = mu / gamma_shape)
  fit <- call_fit(lgspline, fit_formula, dat,
                  list(family = Gamma(link = "log"),
                       keep_weighted_Lambda = TRUE))
  extract_interaction_coverage(fit, "gamma_log", sim_id)
}

scenario_gaussian_ar1 <- function(sim_id) {
  dat <- make_covariates(N, block_size = 5)
  eta <- eta_common(dat, intercept = 0.2, scale = 1)
  err <- draw_block_normal(dat$id, dat$spacetime, block_cor_ar1)
  dat$y <- eta + sigma_gaussian * err
  fit <- call_fit(
    lgspline, fit_formula, dat,
    list(correlation_id = dat$id,
         spacetime = dat$spacetime,
         correlation_structure = "ar(1)",
         VhalfInv_par_init = true_working$ar1)
  )
  rbind(
    extract_interaction_coverage(fit, "gaussian_ar1", sim_id),
    extract_correlation_coverage(fit, "gaussian_ar1", sim_id,
                                 true_working$ar1)
  )
}

scenario_gaussian_matern <- function(sim_id) {
  dat <- make_covariates(N, block_size = 5)
  eta <- eta_common(dat, intercept = 0.2, scale = 1)
  err <- draw_block_normal(dat$id, dat$spacetime, block_cor_matern)
  dat$y <- eta + sigma_gaussian * err
  fit <- call_fit(
    lgspline, fit_formula, dat,
    list(correlation_id = dat$id,
         spacetime = dat$spacetime,
         correlation_structure = "matern",
         VhalfInv_par_init = true_working$matern)
  )
  rbind(
    extract_interaction_coverage(fit, "gaussian_matern", sim_id),
    extract_correlation_coverage(fit, "gaussian_matern", sim_id,
                                 true_working$matern)
  )
}

scenario_poisson_gee_spatial <- function(sim_id) {
  dat <- make_covariates(N, block_size = 5)
  eta <- eta_common(dat, intercept = 1.60, scale = 0.50)
  mu <- exp(eta)
  dat$y <- draw_correlated_count(
    mu, dat$id, dat$spacetime, block_cor_spatial_exponential)
  vf <- make_spatial_exponential_vhalf_functions(dat$id, dat$spacetime)
  fit <- call_fit(
    lgspline, fit_formula, dat,
    list(family = quasipoisson(),
         VhalfInv_fxn = vf$VhalfInv_fxn,
         Vhalf_fxn = vf$Vhalf_fxn,
         VhalfInv_logdet = vf$VhalfInv_logdet,
         VhalfInv_par_init = true_working$spatial_exponential)
  )
  rbind(
    extract_interaction_coverage(fit, "poisson_gee_spatial", sim_id),
    extract_correlation_coverage(fit, "poisson_gee_spatial", sim_id,
                                 true_working$spatial_exponential)
  )
}

scenario_poisson_gee_exchangeable <- function(sim_id) {
  dat <- make_covariates(N, block_size = 5)
  eta <- eta_common(dat, intercept = 1.60, scale = 0.50)
  mu <- exp(eta)
  dat$y <- draw_correlated_count(
    mu, dat$id, dat$spacetime, block_cor_exchangeable)
  vf <- make_exchangeable_vhalf_functions(dat$id)
  fit <- call_fit(
    lgspline, fit_formula, dat,
    list(family = quasipoisson(),
         VhalfInv_fxn = vf$VhalfInv_fxn,
         Vhalf_fxn = vf$Vhalf_fxn,
         VhalfInv_logdet = vf$VhalfInv_logdet,
         VhalfInv_par_init = true_working$exchangeable)
  )
  rbind(
    extract_interaction_coverage(fit, "poisson_gee_exchangeable", sim_id),
    extract_correlation_coverage(fit, "poisson_gee_exchangeable", sim_id,
                                 true_working$exchangeable)
  )
}

scenario_negbin <- function(sim_id) {
  dat <- make_covariates(N)
  eta <- eta_common(dat, intercept = 0.55, scale = 0.65)
  dat$y <- rnbinom(N, size = negbin_size, mu = exp(eta))
  fit <- call_fit(lgspline_negbin, fit_formula, dat)
  extract_interaction_coverage(fit, "negative_binomial", sim_id)
}

scenario_cox <- function(sim_id) {
  dat <- make_covariates(N)
  eta <- eta_common(dat, intercept = 0, scale = 0.70)
  event_time <- rexp(N, rate = 0.12 * exp(eta))
  obs <- apply_censoring(event_time, rate = 0.08)
  dat$time <- obs$y
  fit <- call_fit(lgspline_cox, time_formula, dat,
                  list(status = obs$status))
  extract_interaction_coverage(fit, "cox_ph", sim_id)
}

scenario_weibull_aft <- function(sim_id) {
  dat <- make_covariates(N)
  eta <- eta_common(dat, intercept = 1.5, scale = 0.65)
  event_time <- simulate_weibull_aft(exp(eta), weibull_scale_true)
  obs <- apply_censoring(event_time, rate = 0.06)
  dat$time <- obs$y
  fit <- call_fit(lgspline_weibull, time_formula, dat,
                  list(status = obs$status))
  extract_interaction_coverage(fit, "weibull_aft", sim_id)
}

scenarios <- list(
  gaussian = scenario_gaussian,
  logistic = scenario_logistic,
  poisson = scenario_poisson,
  gamma_log = scenario_gamma_log,
  gaussian_ar1 = scenario_gaussian_ar1,
  gaussian_matern = scenario_gaussian_matern,
  poisson_gee_spatial = scenario_poisson_gee_spatial,
  poisson_gee_exchangeable = scenario_poisson_gee_exchangeable,
  negative_binomial = scenario_negbin,
  cox_ph = scenario_cox,
  weibull_aft = scenario_weibull_aft
)

if(nzchar(scenario_filter)) {
  keep <- trimws(strsplit(scenario_filter, ",")[[1]])
  missing <- setdiff(keep, names(scenarios))
  if(length(missing) > 0) {
    stop("Unknown scenario(s): ", paste(missing, collapse = ", "))
  }
  scenarios <- scenarios[keep]
}

################################################################################

## ## Run Simulation ## ##

cat("lgspline coverage operating-characteristic simulation\n")
cat("N:", N, "\n")
cat("Replicates:", n_sims, "\n")
cat("Fast smoke mode:", fast_smoke, "\n")
if(fast_smoke) cat("Fast smoke K:", fast_K, "\n")
cat("Scenarios:", paste(names(scenarios), collapse = ", "), "\n")
cat("Output CSV:", output_csv, "\n\n")
flush.console()

all_results <- list()
counter <- 0L

start_time <- Sys.time()
for(scenario_name in names(scenarios)) {
  for(sim_id in seq_len(n_sims)) {
    counter <- counter + 1L
    cat(sprintf("[%s] replicate %d/%d\n", scenario_name, sim_id, n_sims))
    flush.console()
    one <- tryCatch(
      scenarios[[scenario_name]](sim_id),
      error = function(e) record_failure(scenario_name, sim_id, e)
    )
    all_results[[counter]] <- one
    write_results(do.call(rbind, all_results), output_csv)
  }
}

results <- do.call(rbind, all_results)
summary_results <- summarize_results(results)

cat("\nRaw result status counts:\n")
print(table(results$scenario, results$status), quote = FALSE)

cat("\nCoverage summary:\n")
print(summary_results, row.names = FALSE)

summary_csv <- write_results(results, output_csv)

cat("\nSaved raw results to:\n", output_csv, "\n", sep = "")
cat("Saved summary results to:\n", summary_csv, "\n", sep = "")
cat("Elapsed time:", format(Sys.time() - start_time), "\n")

invisible(list(raw = results, summary = summary_results))
