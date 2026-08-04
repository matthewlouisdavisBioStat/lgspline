Sys.sleep(5)
strt <- Sys.time()
set.seed(1234)

library(mgcv)
library(dplyr)
library(lgspline)
library(parallel)

## ## Settings ## ##
dims <- c(1, 2, 4)
Ns <- c(500, 5000, 50000)
M_reps <- 24
sigmas <- c(0.5, 2.0)
sp_fixed <- 1e-4
time_cap <- 14400 # time cap in seconds, 4 hours
n_cores <- 12
base_seed <- 1234
tprs_max_knots <- 2000

## ## Conventions ## ##
#  k_per in mgcv is matched to K+1 in lgspline at each q.
#  "Untuned" skips the smoothing-parameter tuning loop. Positive sp in mgcv
#  fixes all smoothing parameters. Each method uses its own minimal ridge.
#  mgcv sp count obtained via gam(..., fit = FALSE), outside the timer.
#  Cubic RS uses bs="cr": cubic regression spline with the integrated-squared
#  -second-derivative ("smoothing spline") penalty, additive + te() for q >= 2.
#  TPRS is a single isotropic thin plate regression spline of all q
#  predictors via s(t1, ..., tq, bs="tp").

################################################################################

## ## Data-Generating Processes ## ##
#  q = 1: Donoho-Johnstone Doppler (1994) in canonical form on [0, 1]:
#         f(t) = 5*sqrt(t*(1-t)) * sin(2*pi*1.05/(t+0.05)).
#  q = 2 and q = 4 use mildly nonstandard predictor distributions: skewed
#  margins, different ranges, weak mixture structure, and correlated predictors
#  for q = 4. The surfaces include localized features and light kinks, but no
#  hard discontinuities.
u_to_t2 <- function(u1, u2) {
  t1 <- 36 * u1^1.65
  t2 <- -3 + 12 * u2^0.62
  cbind(t1, t2)
}

u_to_t4 <- function(u1, u2, u3, u4) {
  t1 <- 55 * u1^1.55
  t2 <- -4 + 16 * u2^0.70
  t3 <- 7 * u3^1.25
  t4 <- -2.2 + 7.0 * u4^0.82 + 0.35 * pmax(u4 - 0.58, 0)^2
  cbind(t1, t2, t3, t4)
}

generate_u2 <- function(n) {
  u1 <- rbeta(n, 1.4, 2.4)
  keep <- runif(n) < 0.25
  u1[keep] <- rbeta(sum(keep), 4.5, 1.8)

  u2 <- rbeta(n, 0.95, 2.2)
  keep <- runif(n) < 0.22
  u2[keep] <- rbeta(sum(keep), 5.5, 2.0)

  cbind(u1, u2)
}

generate_u4 <- function(n) {
  z1 <- rnorm(n)
  z2 <- 0.55 * z1 + sqrt(1 - 0.55^2) * rnorm(n)
  z3 <- -0.35 * z1 + 0.45 * z2 + sqrt(1 - 0.35^2) * rnorm(n)
  z4 <- 0.40 * z1 - 0.25 * z3 + sqrt(1 - 0.40^2) * rnorm(n)

  u1 <- pnorm(z1)
  u2 <- pnorm(z2)
  u3 <- pnorm(z3)
  u4 <- pnorm(z4)

  keep <- runif(n) < 0.20
  u2[keep] <- rbeta(sum(keep), 1.0, 3.2)
  keep <- runif(n) < 0.22
  u4[keep] <- rbeta(sum(keep), 6.0, 2.2)

  cbind(u1, u2, u3, u4)
}

f_q2 <- function(u1, u2) {
  0.75 * sin(2 * pi * u1) * cos(1.5 * pi * u2) +
    1.30 * exp(-55 * (u2 - 0.58 - 0.16 * sin(2.4 * pi * u1))^2) +
    1.35 * exp(-38 * ((u1 - 0.23)^2 + 0.50 * (u2 - 0.78)^2)) -
    1.20 * exp(-70 * ((u1 - 0.72)^2 + (u2 - 0.28)^2)) +
    0.75 * (u1 - 0.5) * (u2 - 0.5) +
    0.22 * abs(u1 - 0.62) -
    0.22 * pmax(u2 - 0.72, 0) +
    0.28 * pmax(u1 + 0.55 * u2 - 1.08, 0)
}

f_q4 <- function(u1, u2, u3, u4) {
  z1 <- u1 - 0.5
  z2 <- u2 - 0.5
  z3 <- u3 - 0.5
  z4 <- u4 - 0.5

  0.65 * sin(2 * pi * u1) +
    0.55 * cos(2 * pi * u2) +
    0.45 * sin(2 * pi * (u3 + 0.25 * u4)) +
    0.65 * z2^2 - 0.45 * z4^2 +
    1.15 * z1 * z2 -
    0.95 * z2 * z3 +
    1.70 * z1 * z3 * z4 +
    1.25 * exp(-55 * ((u1 - 0.75)^2 + (u4 - 0.25)^2)) -
    1.10 * exp(-45 * ((u2 - 0.25)^2 + (u3 - 0.70)^2)) +
    0.18 * abs(u3 - 0.55) +
    0.28 * pmax(u1 - 0.72, 0) * pmax(0.35 - u2, 0) +
    0.16 * pmax(u3 + u4 - 1.15, 0)
}

generate_data <- function(q, n, sigma) {
  if (q == 1) {
    t <- runif(n, 0, 1)
    f_true <- 5*sqrt(t * (1 - t)) * sin(2 * pi * 1.05 / (t + 0.05))
    T_mat <- matrix(t, ncol = 1)
    colnames(T_mat) <- "t1"
  } else if (q == 2) {
    u <- generate_u2(n)
    T_mat <- u_to_t2(u[, 1], u[, 2])
    f_true <- f_q2(u[, 1], u[, 2])
    colnames(T_mat) <- paste0("t", 1:2)
  } else {
    u <- generate_u4(n)
    T_mat <- u_to_t4(u[, 1], u[, 2], u[, 3], u[, 4])
    f_true <- f_q4(u[, 1], u[, 2], u[, 3], u[, 4])
    colnames(T_mat) <- paste0("t", 1:4)
  }
  y <- f_true + rnorm(n, sd = sigma)
  list(T_mat = T_mat, y = y, f_true = f_true)
}

standardization_params <- function(T_mat) {
  ctr <- colMeans(T_mat)
  scl <- apply(T_mat, 2, sd)
  scl[scl == 0] <- 1
  list(center = ctr, scale = scl)
}

apply_standardization <- function(T_mat, params) {
  T_std <- scale(T_mat, center = params$center, scale = params$scale)
  colnames(T_std) <- colnames(T_mat)
  T_std
}

standardize_predictors <- function(T_mat) {
  apply_standardization(T_mat, standardization_params(T_mat))
}

## Knot counts by predictor dimension (used for cubic RS and lgspline)
knots_for <- function(q) switch(as.character(q),
                                "1" = 16, "2" = 8, "4" = 4)

## TPRS basis dimension: mgcv defaults k.def = 8, 27, 100 for d = 1, 2, >=3.
#  For q = 1 we use k_per for parity with the other 1D smoothers.
tprs_k <- function(q, k_per) {
  if (q == 1) k_per
  else if (q == 2) 27
  else 100
}

################################################################################

## ## Fitting Helpers ## ##
## Timed fit with wall-clock cap, returns NULL on error or timeout
timed_fit <- function(fn, cap = time_cap) {
  setTimeLimit(elapsed = cap, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = FALSE))
  tm <- system.time(result <- tryCatch(fn(), error = function(e) {
    message("  error: ", conditionMessage(e))
    NULL
  }))
  setTimeLimit(elapsed = Inf, transient = FALSE)
  if (tm["elapsed"] > cap) result <- NULL
  list(result = result, elapsed = tm["elapsed"])
}

## ## lgspline Fit ## ##
fit_lgspline <- function(T_train, y_train, f_true, q, tuned, k_per) {
  run <- timed_fit(function() {
    lgspline(predictors = T_train, y = y_train, K = k_per * q, opt = tuned,
             include_quartic_terms = TRUE,
             include_2way_interactions = (q >= 2),
             include_quadratic_interactions = (q >= 2),
             include_3way_interactions = (q >= 3),
             return_varcovmat = FALSE, return_G = FALSE,
             return_Ghalf = FALSE, return_U = FALSE,
             estimate_dispersion = FALSE, verbose = FALSE,
             include_warnings = FALSE)
  })
  fit <- run$result
  if (is.null(fit)) return(list(rmse = NA, time = run$elapsed))
  yhat <- tryCatch(predict(fit), error = function(e) rep(NA, length(y_train)))
  list(rmse = sqrt(mean((yhat - f_true)^2, na.rm = TRUE)),
       time = run$elapsed)
}



## ## Cubic Regression Spline
build_formula_cr <- function(q, cnames, k_per) {
  if (q == 1)
    return(as.formula(sprintf("y ~ s(%s, bs = 'cr', k = %d)",
                              cnames[1], k_per)))
  terms <- sprintf("s(%s, bs = 'cr', k = %d)", cnames, k_per)
  k2 <- max(4, floor(sqrt(k_per)))
  for (pr in combn(cnames, 2, simplify = FALSE))
    terms <- c(terms, sprintf("te(%s, %s, bs = 'cr', k = %d)",
                              pr[1], pr[2], k2))
  if (q >= 3) {
    k3 <- max(3, floor(k_per^(1 / 3)))
    for (tr in combn(cnames, 3, simplify = FALSE))
      terms <- c(terms, sprintf("te(%s, %s, %s, bs = 'cr', k = %d)",
                                tr[1], tr[2], tr[3], k3))
  }
  as.formula(paste("y ~", paste(terms, collapse = " + ")))
}

fit_cr <- function(T_train, y_train, f_true, q, tuned, k_per) {
  df <- as.data.frame(T_train)
  df$y <- y_train
  fml <- build_formula_cr(q, colnames(T_train), k_per)

  n_sp <- length(gam(fml, data = df, fit = FALSE)$sp)

  run <- timed_fit(function() {
    if (tuned) {
      gam(fml, data = df, method = "REML")
    } else {
      gam(fml, data = df, method = "REML", sp = rep(sp_fixed, n_sp))
    }
  })
  fit <- run$result
  if (is.null(fit)) return(list(rmse = NA, time = run$elapsed))
  yhat <- predict(fit)
  list(rmse = sqrt(mean((yhat - f_true)^2)), time = run$elapsed)
}



## ## Thin Plate Regression Spline
build_formula_tp <- function(q, cnames, k) {
  vars <- paste(cnames, collapse = ", ")
  as.formula(sprintf("y ~ s(%s, bs = 'tp', k = %d)", vars, k))
}

fit_tp <- function(T_train, y_train, f_true, q, tuned, k_per) {
  df <- as.data.frame(T_train)
  colnames(df) <- colnames(T_train)
  df$y <- y_train

  k <- tprs_k(q, k_per)
  fml <- build_formula_tp(q, colnames(T_train), k)

  ## Pass max.knots via xt to cap basis-construction cost for large N.
  #  Fixed seed inside xt makes the subsample deterministic per call.
  s_args <- list(xt = list(max.knots = tprs_max_knots, seed = 1L))

  ## Count smoothing parameters (one for tp)
  n_sp <- length(gam(fml, data = df, fit = FALSE,
                     xt = s_args$xt)$sp)

  run <- timed_fit(function() {
    if (tuned) {
      gam(fml, data = df, method = "REML",
          xt = s_args$xt)
    } else {
      gam(fml, data = df, method = "REML",
          sp = rep(sp_fixed, n_sp),
          xt = s_args$xt)
    }
  })
  fit <- run$result
  if (is.null(fit)) return(list(rmse = NA, time = run$elapsed))
  yhat <- predict(fit)
  list(rmse = sqrt(mean((yhat - f_true)^2)), time = run$elapsed)
}



## ## Truncated Power Basis ## ##
make_trunc_basis <- function(t_vec, knots, degree = 3) {
  B_poly <- outer(t_vec, 0:degree, "^")
  colnames(B_poly) <- paste0("d", 0:degree)
  B_trunc <- matrix(0, length(t_vec), length(knots))
  for (j in seq_along(knots))
    B_trunc[, j] <- pmax(t_vec - knots[j], 0)^degree
  colnames(B_trunc) <- paste0("tr", seq_along(knots))
  list(poly = B_poly, trunc = B_trunc, full = cbind(B_poly, B_trunc))
}

build_trunc_design <- function(T_mat, K_knots, degree = 3) {
  q <- ncol(T_mat); n <- nrow(T_mat)
  cnames <- colnames(T_mat)

  marginals <- list()
  for (j in seq_len(q)) {
    kn <- quantile(T_mat[, j],
                   probs = seq(0, 1, length.out = K_knots + 2)[-c(1, K_knots + 2)])
    marginals[[j]] <- make_trunc_basis(T_mat[, j], kn, degree)
  }

  X_list <- list(intercept = rep(1, n))
  pen_idx <- c()
  col_count <- 1

  for (j in seq_len(q)) {
    mat_j <- marginals[[j]]$full[, -1, drop = FALSE]
    colnames(mat_j) <- paste0(cnames[j], "_", colnames(marginals[[j]]$full)[-1])
    n_new <- ncol(mat_j)
    new_idx <- col_count + seq_len(n_new)
    X_list <- c(X_list,
                setNames(lapply(seq_len(n_new), function(c) mat_j[, c]),
                         colnames(mat_j)))
    pen_idx <- c(pen_idx, new_idx[(degree + 1):n_new])
    col_count <- col_count + n_new
  }

  if (q >= 2) {
    for (pr in combn(seq_len(q), 2, simplify = FALSE)) {
      j1 <- pr[1]; j2 <- pr[2]
      X_list[[paste0(cnames[j1], "_d1.", cnames[j2], "_d1")]] <-
        marginals[[j1]]$poly[, "d1"] * marginals[[j2]]$poly[, "d1"]
      X_list[[paste0(cnames[j1], "_d2.", cnames[j2], "_d1")]] <-
        marginals[[j1]]$poly[, "d2"] * marginals[[j2]]$poly[, "d1"]
      X_list[[paste0(cnames[j1], "_d1.", cnames[j2], "_d2")]] <-
        marginals[[j1]]$poly[, "d1"] * marginals[[j2]]$poly[, "d2"]
      col_count <- col_count + 3
    }
  }
  if (q >= 3) {
    for (tr in combn(seq_len(q), 3, simplify = FALSE)) {
      X_list[[paste0(cnames[tr[1]], "_d1.", cnames[tr[2]], "_d1.",
                     cnames[tr[3]], "_d1")]] <-
        marginals[[tr[1]]]$poly[, "d1"] *
        marginals[[tr[2]]]$poly[, "d1"] *
        marginals[[tr[3]]]$poly[, "d1"]
      col_count <- col_count + 1
    }
  }

  X <- do.call(cbind, X_list)
  colnames(X) <- names(X_list)
  list(X = X, pen_idx = pen_idx)
}

fit_trunc_poly <- function(T_train, y_train, f_true, q, tuned, k_per) {
  des <- build_trunc_design(T_train, k_per, degree = 4)
  X <- des$X; p <- ncol(X)

  run <- timed_fit(function() {
    if (tuned) {
      S_mat <- matrix(0, p, p)
      diag(S_mat)[des$pen_idx] <- 1
      fit <- magic(y_train, X, sp = -1, S = list(S_mat), off = 1)
      fit$b
    } else {
      pen <- rep(0, p)
      pen[des$pen_idx] <- sp_fixed
      XtX <- crossprod(X)
      diag(XtX) <- diag(XtX) + pen
      as.numeric(solve(XtX, crossprod(X, y_train)))
    }
  })
  coefs <- run$result
  if (is.null(coefs)) return(list(rmse = NA, time = run$elapsed))
  yhat <- X %*% coefs
  list(rmse = sqrt(mean((yhat - f_true)^2)), time = run$elapsed)
}

################################################################################

## ## Parallel Setup ## ##
scenarios <- expand.grid(q = dims, N = Ns, sigma = sigmas,
                         stringsAsFactors = FALSE)
scenarios$K <- sapply(scenarios$q, knots_for)

run_one_rep <- function(args) {
  q <- args$q
  N <- args$N
  sigma <- args$sigma
  k_per <- args$K
  b <- args$rep

  suppressMessages({
    library(mgcv)
    library(lgspline)
  })

  train <- generate_data(q, N, sigma)
  T_fit <- standardize_predictors(train$T_mat)
  rows <- list()

  for (tuned in c(TRUE, FALSE)) {
    tune_lab <- ifelse(tuned, "tuned", "untuned")

    res_lg <- tryCatch(
      fit_lgspline(T_fit, train$y, train$f_true, q, tuned, k_per),
      error = function(e) list(rmse = NA, time = NA))

    res_cr <- tryCatch(
      fit_cr(T_fit, train$y, train$f_true, q, tuned, k_per),
      error = function(e) list(rmse = NA, time = NA))

    res_tp <- tryCatch(
      fit_tp(T_fit, train$y, train$f_true, q, tuned, k_per),
      error = function(e) list(rmse = NA, time = NA))

    res_tr <- tryCatch(
      fit_trunc_poly(T_fit, train$y, train$f_true, q, tuned, k_per),
      error = function(e) list(rmse = NA, time = NA))

    rows[[length(rows) + 1]] <- data.frame(
      q = q, N = N, sigma = sigma, K = k_per,
      tuning = tune_lab, rep = b,
      method = c("lgspline", "Cubic RS", "TPRS", "Trunc. poly"),
      rmse = c(res_lg$rmse, res_cr$rmse, res_tp$rmse, res_tr$rmse),
      time = c(res_lg$time, res_cr$time, res_tp$time, res_tr$time),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

## ## Tasks and Execution ## ##
tasks <- list()
for (i in seq_len(nrow(scenarios))) {
  for (m in seq_len(M_reps)) {
    tasks[[length(tasks) + 1]] <- c(as.list(scenarios[i, ]), list(rep = m))
  }
}

task_seeds <- base_seed + seq_along(tasks)
run_seeded <- function(idx) {
  set.seed(task_seeds[idx])
  run_one_rep(tasks[[idx]])
}

cl <- makeCluster(n_cores)
clusterExport(cl, c("tasks", "task_seeds", "run_one_rep",
                    "generate_data", "knots_for", "tprs_k", "tprs_max_knots",
                    "u_to_t2", "u_to_t4", "generate_u2", "generate_u4",
                    "f_q2", "f_q4",
                    "standardization_params", "apply_standardization",
                    "standardize_predictors",
                    "timed_fit", "time_cap", "sp_fixed",
                    "fit_lgspline", "fit_cr", "fit_tp", "fit_trunc_poly",
                    "build_formula_cr", "build_formula_tp",
                    "make_trunc_basis", "build_trunc_design"))

results_list <- parLapply(cl, seq_along(tasks), run_seeded)
stopCluster(cl)

results <- do.call(rbind, results_list)
session_info <- sessionInfo()
print(session_info)
save(results, session_info, file = "performance_results.RData")

################################################################################

## ## Summarize Results

load('performance_results.RData')
## ## Summaries ## ##
summary_tab <- results %>%
  group_by(q, N, sigma, K, tuning, method) %>%
  summarise(
    rmse_mean = mean(rmse, na.rm = TRUE),
    rmse_sd = sd(rmse, na.rm = TRUE),
    # Timing over completed only (methods reached timing cap)
    time_mean_completed = mean(time[!is.na(rmse)], na.rm = TRUE),
    time_sd_completed   = sd(time[!is.na(rmse)], na.rm = TRUE),
    # Timing over all
    time_mean_all = mean(pmin(time, time_cap), na.rm = TRUE),
    n_total = n(),
    n_ok = sum(!is.na(rmse)),
    n_timeout = sum(is.na(rmse) & !is.na(time)),
    .groups = "drop")

cat("\n\nSummary (mean +/- SD):\n")
print(as.data.frame(summary_tab), digits = 3, max = 999999)
cat(sprintf("\nTotal elapsed: %.1f min\n",
            as.numeric(difftime(Sys.time(), strt, units = "mins"))))

## ## Faceted Timing Figure ## ##
#  Log-log time vs. N, faceted by (q, sigma), colored by method, tuned vs. untuned
#  Cap-exceeded cells are imputed at the cap value and marked with a triangle.
#  All other cells plotted as-is; accuracy caveats handled in text.

## Aggregate to cell means and SDs
library(dplyr)
plot_tab <- results %>%
  group_by(q, N, sigma, tuning, method) %>%
  summarise(
    time_mean = mean(time, na.rm = TRUE),
    time_sd   = sd(time, na.rm = TRUE),
    n_ok      = sum(!is.na(rmse)),
    .groups = 'drop'
  ) %>%
  mutate(is_cap = (n_ok < 12 & time_mean > 0.9 * time_cap))

## Method colors and tuning linetypes
methods <- c('lgspline', 'Cubic RS', 'TPRS', 'Trunc. poly')
cols <- c(lgspline      = 'black',
          `Cubic RS`    = 'red',
          TPRS          = 'blue',
          `Trunc. poly` = 'darkgreen')
ltys <- c(tuned = 1, untuned = 2)

png("timing_facets.png",
    width = 1800,
    height = 1200,
    res = 150)

## Layout: 2x3 panels on top, legend strip on bottom
layout(matrix(c(1, 2, 3,
                4, 5, 6,
                7, 7, 7),
              nrow = 3, byrow = TRUE),
       heights = c(1, 1, 0.22))
par(mar = c(4.5, 4.8, 3, 1), oma = c(0, 0, 3, 0))

## Log-scale y-axis ticks (powers of 10)
y_major <- 10^seq(-3, 4, by = 1)
y_major_labs <- c(
  expression(10^-3), expression(10^-2), expression(10^-1),
  expression(10^0), expression(10^1), expression(10^2),
  expression(10^3), expression(10^4)
)
y_minor <- as.vector(outer(2:9, 10^seq(-3, 3, by = 1)))

for (sig in c(0.5, 2.0)) {
  for (qv in c(1, 2, 4)) {
    sub <- plot_tab[plot_tab$q == qv & plot_tab$sigma == sig, ]

    plot(NA, NA,
         log = 'xy',
         xlim = range(plot_tab$N),
         ylim = c(0.001, 10000),
         xlab = 'N',
         ylab = 'Time (s, log scale)',
         main = bquote(q == .(qv) * ',' ~ sigma == .(sig)),
         axes = FALSE,
         cex.lab = 1.2, cex.main = 1.3)

    abline(h = y_major, col = 'grey88', lwd = 0.7)
    abline(h = y_minor, col = 'grey96', lwd = 0.4)

    axis(1, cex.axis = 1.0)
    axis(2, at = y_major, labels = y_major_labs, las = 1, cex.axis = 1.0)
    axis(2, at = y_minor, labels = FALSE, tcl = -0.25)
    box()

    for (m in methods) {
      for (tn in c('tuned', 'untuned')) {
        dd <- sub[sub$method == m & sub$tuning == tn, ]
        dd <- dd[order(dd$N), ]
        if (nrow(dd) == 0) next

        ## Impute cap-exceeded at the cap value
        y_plot <- ifelse(dd$is_cap, time_cap, dd$time_mean)

        lines(dd$N, y_plot, col = cols[m], lty = ltys[tn], lwd = 2)

        ## Normal filled circles
        ok <- !dd$is_cap
        points(dd$N[ok], y_plot[ok], col = cols[m], pch = 19, cex = 1.1)

        ## Cap-exceeded triangles
        points(dd$N[dd$is_cap], y_plot[dd$is_cap],
               col = cols[m], pch = 17, cex = 1.6)
      }
    }
  }
}

## Outer title
mtext('Wall-clock time by sample size, method, and tuning regime',
      outer = TRUE, cex = 1.3, font = 2)

## Legend strip: two legends with a gap between them
par(mar = c(0, 0, 0, 0))
plot.new()

legend(
  x = 0.28, y = 0.5,
  legend = methods,
  col = cols[methods],
  lty = 1,
  pch = 19,
  lwd = 2,
  bty = 'n',
  cex = 1.4,
  horiz = TRUE,
  xjust = 0.5, yjust = 0.5
)

legend(
  x = 0.78, y = 0.5,
  legend = c('tuned', 'untuned'),
  col = c('grey30', 'grey30'),
  lty = c(1, 2),
  lwd = 2,
  bty = 'n',
  cex = 1.4,
  horiz = TRUE,
  xjust = 0.5, yjust = 0.5
)

dev.off()

################################################################################

## ## Doppler Fit Comparison Figure ## ##
#  One replicate at N=5000, sigma=0.5, tuned. Overlay truth + 4 methods.

set.seed(1234)
dop <- generate_data(q = 1, n = 5000, sigma = 0.5)
t_vec <- as.numeric(dop$T_mat[, 1])
y_vec <- dop$y
dop_std <- standardization_params(dop$T_mat)
T_dop_fit <- apply_standardization(dop$T_mat, dop_std)

## Dense grid for plotting true + predicted curves
t_grid <- seq(0, 1, length.out = 2000)
t_grid_mat <- matrix(t_grid, ncol = 1, dimnames = list(NULL, 't1'))
t_grid_fit <- apply_standardization(t_grid_mat, dop_std)
f_grid <- 5*sqrt(t_grid * (1 - t_grid)) * sin(2 * pi * 1.05 / (t_grid + 0.05))

## lgspline
fit_lg <- lgspline(predictors = T_dop_fit, y = y_vec, K = 16,
                   include_quartic_terms = TRUE,
                   return_varcovmat = FALSE, return_G = FALSE,
                   return_Ghalf = FALSE, return_U = FALSE,
                   estimate_dispersion = FALSE,
                   verbose = FALSE, include_warnings = FALSE)
pred_lg <- as.numeric(predict(fit_lg, newdata = t_grid_fit))

## Cubic RS (mgcv, smoothing spline penalty via bs="cr")
df_sim <- data.frame(t1 = as.numeric(T_dop_fit[, 1]), y = y_vec)
fit_cr_dop <- gam(y ~ s(t1, bs = 'cr', k = 16), data = df_sim, method = 'REML')
pred_cr <- predict(fit_cr_dop, newdata = data.frame(t1 = as.numeric(t_grid_fit[, 1])))

## TPRS (mgcv, isotropic thin plate via bs="tp"); 1D so just a single smoother
fit_tp_dop <- gam(y ~ s(t1, bs = 'tp', k = 16), data = df_sim, method = 'REML')
pred_tp <- predict(fit_tp_dop, newdata = data.frame(t1 = as.numeric(t_grid_fit[, 1])))

## Truncated polynomial: stack training + grid so knots are shared
t_stack <- rbind(T_dop_fit, t_grid_fit)
des_stack <- build_trunc_design(t_stack, K_knots = 16, degree = 4)
X_stack <- des_stack$X
p <- ncol(X_stack)

N_train <- nrow(dop$T_mat)
X_train <- X_stack[seq_len(N_train), , drop = FALSE]
X_grid  <- X_stack[(N_train + 1):nrow(X_stack), , drop = FALSE]

S_mat <- matrix(0, p, p)
diag(S_mat)[des_stack$pen_idx] <- 1
fit_trp <- magic(y_vec, X_train, sp = -1, S = list(S_mat), off = 1)
pred_trp <- as.numeric(X_grid %*% fit_trp$b)

## Plot
png("doppler_fits.png",
    width = 1600,
    height = 900,
    res = 150)

par(mar = c(5, 5, 3, 2))
plot(t_grid, f_grid,
     type = 'l',
     col = 'gold',
     lwd = 3,
     xlab = 't',
     ylab = 'f(t)',
     main = 'Doppler fits (N = 5000, sigma = 0.5)',
     cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.4,
     ylim = range(c(f_grid, pred_lg, pred_cr, pred_tp, pred_trp),
                  na.rm = TRUE))

## Scatter of noisy data, subsampled for legibility
pts_idx <- seq_along(t_vec)
points(t_vec[pts_idx], y_vec[pts_idx],
       col = rgb(0.5, 0.5, 0.5, 0.3), pch = 19, cex = 0.5)

lines(t_grid, pred_lg,  col = 'black',     lwd = 2, lty = 1)
lines(t_grid, pred_cr,  col = 'red',       lwd = 2, lty = 2)
lines(t_grid, pred_tp,  col = 'blue',      lwd = 2, lty = 3)
lines(t_grid, pred_trp, col = 'darkgreen', lwd = 2, lty = 4)

legend('topright',
       legend = c('truth', 'lgspline', 'Cubic RS', 'TPRS', 'Trunc. poly'),
       col = c('gold', 'black', 'red', 'blue', 'darkgreen'),
       lty = c(1, 1, 2, 3, 4),
       lwd = c(3, 2, 2, 2, 2),
       bty = 'n', cex = 1.1)

dev.off()
