library(lgspline)
library(splines)
set.seed(1234)

## True surface defined in centered normalized coordinates.
#  t10, t20 are the domain lower corners; width is the side length.
f_true <- function(t1, t2, t10, t20, width) {
  t1r <- (t1 - (t10 + width / 2)) / width
  t2r <- (t2 - (t20 + width / 2)) / width
  2 + 0.8 * sin(4 * t1r) + 0.7 * cos(4 * t2r) +
    0.6 * t1r * t2r - 0.4 * t1r^2 + 0.3 * t2r^2 +
    0.15 * t1r^3 - 0.10 * t2r^3
}

## Extract the mean condition number across partitions from a fitted lgspline.
#  model$G is a list of per-partition penalized Gram matrices.
mean_kappa <- function(model) {
  kappas <- unlist(lapply(model$G, kappa))
  mean(kappas, na.rm = TRUE)
}

## Build truncated power basis (cubic with one knot)
tp1d <- function(x, knot) {
  cbind(1, x, x^2, x^3, pmax(x - knot, 0)^3)
}

## Run one Monte Carlo replication for a given wiggle_penalty level.
#  Also computes B-spline and truncated power baselines on same data.
run_rep <- function(t10, t20, width, N, noise_sd,
                    wiggle_penalty, flat_ridge_penalty = 1) {
  
  t1_range <- c(t10, t10 + width)
  t2_range <- c(t20, t20 + width)
  
  t1 <- runif(N, t1_range[1], t1_range[2])
  t2 <- runif(N, t2_range[1], t2_range[2])
  Tmat <- cbind(t1, t2)
  
  y <- f_true(t1, t2, t10, t20, width) + rnorm(N, sd = noise_sd)
  
  ## lgspline fit
  fit <- tryCatch(
    lgspline(
      Tmat, y,
      unique_penalty_per_partition = FALSE,
      unique_penalty_per_predictor = FALSE,
      wiggle_penalty = wiggle_penalty,
      include_quadratic_interactions = TRUE,
      include_quartic_terms = FALSE,
      K = 3,
      opt = FALSE
    ),
    error = function(e) NULL
  )
  
  ## B-spline baseline
  Bt1 <- bs(t1, df = 6, degree = 3, intercept = TRUE,
            Boundary.knots = t1_range)
  Bt2 <- bs(t2, df = 6, degree = 3, intercept = TRUE,
            Boundary.knots = t2_range)
  
  X_bs <- matrix(0, nrow = N, ncol = ncol(Bt1) * ncol(Bt2))
  for (i in seq_len(N)) {
    X_bs[i, ] <- as.vector(outer(Bt1[i, ], Bt2[i, ]))
  }
  fit_bs <- lm.fit(X_bs, y)
  
  ## Truncated power baseline
  k1 <- mean(t1_range)
  k2 <- mean(t2_range)
  
  Tt1 <- tp1d(t1, k1)
  Tt2 <- tp1d(t2, k2)
  
  X_tp <- matrix(0, nrow = N, ncol = ncol(Tt1) * ncol(Tt2))
  for (i in seq_len(N)) {
    X_tp[i, ] <- as.vector(outer(Tt1[i, ], Tt2[i, ]))
  }
  fit_tp <- lm.fit(X_tp, y)
  
  if (is.null(fit) ||
      any(!is.finite(fit_bs$coefficients)) ||
      any(!is.finite(fit_tp$coefficients))) {
    return(c(kappa_mean = NA, rmse_truth = NA,
             kappa_bs = NA, rmse_bs = NA,
             kappa_tp = NA, rmse_tp = NA))
  }
  
  ## Grid for truth RMSE
  t1g <- seq(t1_range[1], t1_range[2], length.out = 81)
  t2g <- seq(t2_range[1], t2_range[2], length.out = 81)
  grd <- expand.grid(t1 = t1g, t2 = t2g)
  
  y_truth <- f_true(grd$t1, grd$t2, t10, t20, width)
  
  ## lgspline prediction
  pred_g <- tryCatch(
    predict(fit, newdata = as.matrix(grd)),
    error = function(e) rep(NA_real_, nrow(grd))
  )
  
  ## B-spline prediction
  Bgt1 <- bs(grd$t1, df = 6, degree = 3, intercept = TRUE,
             Boundary.knots = t1_range)
  Bgt2 <- bs(grd$t2, df = 6, degree = 3, intercept = TRUE,
             Boundary.knots = t2_range)
  
  Xg_bs <- matrix(0, nrow = nrow(grd), ncol = ncol(X_bs))
  for (i in seq_len(nrow(grd))) {
    Xg_bs[i, ] <- as.vector(outer(Bgt1[i, ], Bgt2[i, ]))
  }
  pred_bs <- as.vector(Xg_bs %*% fit_bs$coefficients)
  
  ## Truncated power prediction
  Tg1 <- tp1d(grd$t1, k1)
  Tg2 <- tp1d(grd$t2, k2)
  
  Xg_tp <- matrix(0, nrow = nrow(grd), ncol = ncol(X_tp))
  for (i in seq_len(nrow(grd))) {
    Xg_tp[i, ] <- as.vector(outer(Tg1[i, ], Tg2[i, ]))
  }
  pred_tp <- as.vector(Xg_tp %*% fit_tp$coefficients)
  
  if (!all(is.finite(pred_g)) ||
      !all(is.finite(pred_bs)) ||
      !all(is.finite(pred_tp))) {
    return(c(kappa_mean = NA, rmse_truth = NA,
             kappa_bs = NA, rmse_bs = NA,
             kappa_tp = NA, rmse_tp = NA))
  }
  
  c(
    kappa_mean = mean_kappa(fit),
    rmse_truth = sqrt(mean((pred_g - y_truth)^2)),
    
    kappa_bs = kappa(crossprod(X_bs)),
    rmse_bs = sqrt(mean((pred_bs - y_truth)^2)),
    
    kappa_tp = kappa(crossprod(X_tp)),
    rmse_tp = sqrt(mean((pred_tp - y_truth)^2))
  )
}

## Grid of penalties
wiggle_grid <- c(0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1)
n_reps <- 250

mc_results <- lapply(wiggle_grid, function(wp) {
  
  reps <- t(replicate(n_reps, run_rep(
    t10 = 0, t20 = 0, width = 1,
    N = 1000, noise_sd = 1,
    wiggle_penalty = wp
  )))
  
  data.frame(
    wiggle_penalty = wp,
    
    ## lgspline
    kappa_mean = mean(reps[, "kappa_mean"], na.rm = TRUE),
    kappa_sd   = sd(reps[, "kappa_mean"], na.rm = TRUE),
    
    rmse_mean  = mean(reps[, "rmse_truth"], na.rm = TRUE),
    rmse_sd    = sd(reps[, "rmse_truth"], na.rm = TRUE),
    
    rmse_cv = sd(reps[, "rmse_truth"], na.rm = TRUE) /
      mean(reps[, "rmse_truth"], na.rm = TRUE),
    
    ## B-spline
    kappa_bs_mean = mean(reps[, "kappa_bs"], na.rm = TRUE),
    kappa_bs_sd   = sd(reps[, "kappa_bs"], na.rm = TRUE),
    
    rmse_bs_mean  = mean(reps[, "rmse_bs"], na.rm = TRUE),
    rmse_bs_sd    = sd(reps[, "rmse_bs"], na.rm = TRUE),
    
    ## Truncated power
    kappa_tp_mean = mean(reps[, "kappa_tp"], na.rm = TRUE),
    kappa_tp_sd   = sd(reps[, "kappa_tp"], na.rm = TRUE),
    
    rmse_tp_mean  = mean(reps[, "rmse_tp"], na.rm = TRUE),
    rmse_tp_sd    = sd(reps[, "rmse_tp"], na.rm = TRUE)
  )
})

summary_mc <- do.call(rbind, mc_results)
rownames(summary_mc) <- NULL

save(summary_mc, file = "stability_results.RData")
print(signif(summary_mc, 4))

################################################################################

## ## Dual-Axis Figure: Condition Number and RMSE vs. Penalty ## ##
## lgspline kappa (log scale, left) and RMSE (linear scale, right) vs. lambda_w
#  B-spline and truncated power shown as horizontal reference lines

load('stability_results.RData')

## x positions for display: lambda=0 placed at a visible left slot
lam <- summary_mc$wiggle_penalty
x_pos <- ifelse(lam == 0, 1e-12, lam)   # push zero onto log axis at 1e-12

png("conditioning.png",
    width = 1400,
    height = 900,
    res = 150)

par(mar = c(5, 5, 3, 5))

## Left axis: log-scale kappa
plot(x_pos, summary_mc$kappa_mean,
     log = 'xy',
     type = 'b',
     pch = 19,
     lwd = 2,
     col = 'black',
     xlab = expression(lambda[w]),
     ylab = expression('Condition number '~kappa(bold(G)[k])),
     xaxt = 'n',
     ylim = c(1e3, 1e12),
     cex = 1.4,
     cex.lab = 1.3,
     cex.axis = 1.1)

## Custom x-axis with lambda=0 label
axis(1,
     at = c(1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1),
     labels = c('0', expression(10^-10), expression(10^-8),
                expression(10^-6), expression(10^-4),
                expression(10^-2), '1'),
     cex.axis = 1.1)

## SD bands for kappa
kappa_lo <- pmax(summary_mc$kappa_mean - summary_mc$kappa_sd, 1e2)
kappa_hi <- summary_mc$kappa_mean + summary_mc$kappa_sd
segments(x_pos, kappa_lo, x_pos, kappa_hi, col = 'black', lwd = 1)

## Reference lines for B-spline and truncated power
abline(h = mean(summary_mc$kappa_bs_mean), col = 'blue', lty = 3, lwd = 2)
abline(h = mean(summary_mc$kappa_tp_mean), col = 'red',  lty = 3, lwd = 2)
text(1e-11, mean(summary_mc$kappa_bs_mean) * 2.5,
     'B-spline', col = 'blue', pos = 4, cex = 1.1)
text(1e-11, mean(summary_mc$kappa_tp_mean) * 0.4,
     'Trunc. power', col = 'red', pos = 4, cex = 1.1)

## Right axis: RMSE
par(new = TRUE)
rmse_range <- range(summary_mc$rmse_mean)
plot(x_pos, summary_mc$rmse_mean,
     log = 'x',
     type = 'b',
     pch = 17,
     lty = 2,
     lwd = 2,
     col = 'darkgreen',
     axes = FALSE,
     xlab = '',
     ylab = '',
     ylim = c(0.10, 0.23),
     cex = 1.4)
axis(4, cex.axis = 1.1, col = 'darkgreen', col.axis = 'darkgreen')
mtext('RMSE', side = 4, line = 3, col = 'darkgreen', cex = 1.3)

## SD bands for RMSE
rmse_lo <- summary_mc$rmse_mean - summary_mc$rmse_sd
rmse_hi <- summary_mc$rmse_mean + summary_mc$rmse_sd
segments(x_pos, rmse_lo, x_pos, rmse_hi, col = 'darkgreen', lwd = 1)

## Annotate the sweet spot at lambda_w = 1e-2
idx_sweet <- which(lam == 1e-2)
points(x_pos[idx_sweet], summary_mc$rmse_mean[idx_sweet],
       pch = 1, cex = 3, lwd = 2, col = 'darkgreen')

## Legend
legend('topleft',
       legend = c(expression(kappa(bold(G)[k])~' (lgspline)'),
                  'RMSE (lgspline)'),
       col = c('black', 'darkgreen'),
       pch = c(19, 17),
       lty = c(1, 2),
       lwd = 2,
       bty = 'n',
       cex = 1.1)

dev.off()