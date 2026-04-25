set.seed(1234)
library(lgspline)
library(parallel)

## Prep volcano data
data('volcano')
volcano_long <-
  Reduce('rbind', lapply(1:nrow(volcano), function(i){
    t(sapply(1:ncol(volcano), function(j){
      c(i, j, volcano[i, j])
    }))
  }))
colnames(volcano_long) <- c('Length', 'Width', 'Height')

predictors <- volcano_long[, c(1, 2)]
response   <- volcano_long[, 3]
print(dim(volcano_long))
# [1] 5307    3

## Holdout split (10%)
set.seed(1234)
N <- nrow(predictors)
test_idx <- sample.int(N, size = floor(0.1 * N))
train_idx  <- setdiff(seq_len(N), test_idx)

predictors_train <- predictors[train_idx, , drop = FALSE]
response_train <- response[train_idx]

predictors_test <- predictors[test_idx, , drop = FALSE]
response_test <- response[test_idx]

## Partition counts to test
K_vals <- c(50, 100, 250, 500, 1000) - 1
n_cores <- 12
B_reps <- 5

## Common fit arguments (no tuning)
fit_args <- list(
  opt = FALSE,
  include_quadratic_interactions = TRUE,
  include_constrain_first_deriv = FALSE,
  include_constrain_second_deriv = FALSE,
  parallel_eigen = TRUE,
  parallel_find_neighbors = TRUE,
  parallel_make_constraint = TRUE,
  return_varcovmat = FALSE,
  return_G = FALSE,
  return_Ghalf = FALSE,
  return_U = FALSE,
  estimate_dispersion = FALSE
)

results <- data.frame()
results_opt <- data.frame()

safe_fit <- function(expr, max_tries = 3) {
  for (i in 1:max_tries) {
    fit <- tryCatch(eval(expr), error = function(e) NULL)
    if (!is.null(fit)) return(fit)
    cat(sprintf("    try %d failed, retrying...\n", i))
  }
  return(NULL)
}

safe_predict <- function(fit_expr, pred_expr, max_tries = 3) {
  for (i in 1:max_tries) {
    fit <- safe_fit(fit_expr, max_tries = 1)
    if (!is.null(fit)) {
      yhat <- tryCatch(eval(pred_expr), error = function(e) NULL)
      if (!is.null(yhat)) return(list(fit = fit, yhat = yhat))
    }
    cat(sprintf("    predict try %d failed, refitting...\n", i))
  }
  return(list(fit = NULL, yhat = NULL))
}

for (K in K_vals) {
  cat(sprintf("\nK = %d\n", K))
  
  ## Parallel fits
  cl <- makeCluster(n_cores)
  times_par <- numeric(B_reps)
  rmse_par <- numeric(B_reps)
  for (b in 1:B_reps) {
    args <- c(list(predictors = predictors_train, y = response_train, K = K, cl = cl),
              fit_args)
    tm <- system.time({
      res <- safe_predict(
        quote(do.call(lgspline, args)),
        quote(predict(fit, new_predictors = predictors_test))
      )
      fit <- res$fit
      yhat <- res$yhat
    })
    times_par[b] <- tm["elapsed"]
    if (!is.null(yhat)) {
      rmse_par[b] <- sqrt(mean((yhat - response_test)^2))
    } else {
      rmse_par[b] <- NA
    }
    cat(sprintf("  parallel rep %d: %.1fs  RMSE=%.3f\n", b, tm["elapsed"], rmse_par[b]))
  }
  stopCluster(cl)
  
  ## Sequential fits
  times_seq <- numeric(B_reps)
  rmse_seq <- numeric(B_reps)
  for (b in 1:B_reps) {
    args    <- c(list(predictors = predictors_train, y = response_train, K = K), fit_args)
    args$cl <- NULL
    Sys.sleep(0.5)
    tm <- system.time({
      res <- safe_predict(
        quote(do.call(lgspline, args)),
        quote(predict(fit, new_predictors = predictors_test))
      )
      fit <- res$fit
      yhat <- res$yhat
    })
    times_seq[b] <- tm["elapsed"]
    if (!is.null(yhat)) {
      rmse_seq[b] <- sqrt(mean((yhat - response_test)^2))
    } else {
      rmse_seq[b] <- NA
    }
    cat(sprintf("  sequential rep %d: %.1fs  RMSE=%.3f\n", b, tm["elapsed"], rmse_seq[b]))
  }
  
  speedup <- times_seq / times_par
  
  ## Collect results
  results <- rbind(results, data.frame(
    K = K,
    partitions = K + 1,
    time_par_mean = mean(times_par),
    time_par_sd = sd(times_par),
    time_seq_mean = mean(times_seq),
    time_seq_sd = sd(times_seq),
    speedup_mean = mean(speedup),
    speedup_sd = sd(speedup),
    rmse_mean = mean(rmse_par, na.rm = TRUE),
    rmse_sd = sd(rmse_par, na.rm = TRUE)
  ))
}

save(results, file = 'volcano_results.RData')

## Print results
cat("\n\nResults:\n")
print(results, digits = 3)

################################################################################

## ## Speedup Curve Figure ## ##
#  Speedup = sequential time / parallel time, with error bars over B_reps
load('volcano_results.RData')

png("speedup.png",
    width = 1200,
    height = 900,
    res = 150)

par(mar = c(5, 5, 3, 2))

## Propagate SD via the ratio: sd(seq/par) ~= |seq/par| * sqrt((sd_s/mean_s)^2 + (sd_p/mean_p)^2)
cv_seq <- results$time_seq_sd / results$time_seq_mean
cv_par <- results$time_par_sd / results$time_par_mean
speedup_sd <- results$speedup * sqrt(cv_seq^2 + cv_par^2)

plot(results$partitions, results$speedup_mean,
     log = 'x',
     type = 'b',
     pch = 19,
     lwd = 2,
     col = 'black',
     xlab = expression('Partitions ('*K*'+1)'),
     ylab = expression('Speedup ('*t[seq]/t[par]*')'),
     ylim = c(0, max(results$speedup_mean + results$speedup_sd, na.rm = TRUE) * 1.15),
     cex = 1.5,
     cex.lab = 1.3,
     cex.axis = 1.1,
     main = 'Parallel speedup on the volcano dataset (12 cores)',
     cex.main = 1.3)

## Error bars
segments(results$partitions,
         pmax(results$speedup_mean - results$speedup_sd, 0),
         results$partitions,
         results$speedup_mean + results$speedup_sd,
         lwd = 1)

## Break-even line and crossover annotation
abline(h = 1, lty = 2, lwd = 1.5, col = 'grey40')
text(min(results$partitions) * 1.2, 1.08, '',
     col = 'grey40', pos = 4, cex = 1.0)

dev.off()

################################################################################

## Final fit for illustration, K + 1 = 250 with tuning and plotting
set.seed(1234)
cl <- makeCluster(12)
tm_par <-
  system.time({
    final_fit <- lgspline(
        predictors,
        response,
        K = 249,
        cl = cl,
        include_quartic_terms = FALSE,
        include_constrain_first_deriv = FALSE,
        include_constrain_second_deriv = FALSE,
        parallel_eigen = TRUE,
        parallel_find_neighbors = TRUE,
        parallel_make_constraint = TRUE,
        return_varcovmat = FALSE,
        return_G = FALSE,
        return_Ghalf = FALSE,
        return_U = FALSE,
        estimate_dispersion = FALSE,
        initial_wiggle = 1e-8,
        initial_flat = 1e-8
      )})
stopCluster(cl)

## Plotting on new data with interactive visual + formulas
new_input <- expand.grid(seq(min(volcano_long[,1]),
                             max(volcano_long[,1]),
                             length.out = 350),
                         seq(max(volcano_long[,2]),
                             min(volcano_long[,2]),
                             length.out = 350))
plot(final_fit,
     new_predictors = new_input,
     show_formulas = TRUE,
     custom_response_lab = "Height",
     custom_title = 'Volcano 3-D Map',
     text_size_formula = 14,
     digits = 4)
tm_par
# user  system elapsed 
# 51.27    8.74   80.83 