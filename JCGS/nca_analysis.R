library(lgspline)
set.seed(1234)

## Data
data("Theoph")

## Set up base data frame
#  Raw, natural scale
df <- Theoph[,c('Time','Dose','conc','Subject')]
levels(df$Subject) <- as.numeric(unique(as.character(df$Subject)))
df <- df[order(df$Time, df$Dose),]

################################################################################

## ## Correlation Helpers ## ##
## Compute subject ordering and block sizes from the reordered data frame
subject_order <- as.integer(as.character(df$Subject))   # subject label per row
unique_subjects <- unique(subject_order)                # in appearance order
n_blocks <- length(unique_subjects)

## block_sizes[k] = number of rows belonging to the k-th subject (in row order)
block_sizes <- sapply(unique_subjects, function(s) sum(subject_order == s))
N <- nrow(df)

## Build a tridiagonal correlation matrix for a block of size n
tridiag_corr <- function(n, rho) {
  C <- diag(n)
  if (n > 1) {
    for (i in 1:(n - 1)) {
      C[i, i + 1] <- rho
      C[i + 1, i] <- rho
    }
  }
  C
}

## V^{-1/2}: block-diagonal, one block per subject
VhalfInv_fxn <- function(par) {
  rho <- tanh(par[1])
  block_list <- lapply(block_sizes, function(n) {
    matinvsqrt(tridiag_corr(n, rho))
  })
  ## Assemble block-diagonal matrix
  collapse_block_diagonal(block_list)     # sparse; lgspline accepts Matrix objects
}

## V^{1/2}: needed because we are fitting a Gaussian GLM (identity link)
Vhalf_fxn <- function(par) {
  rho <- tanh(par[1])
  block_list <- lapply(block_sizes, function(n) {
    matsqrt(tridiag_corr(n, rho))
  })
  collapse_block_diagonal(block_list)
}

## Log-determinant of V^{-1/2} (sum over blocks for efficiency)
VhalfInv_logdet <- function(par) {
  rho <- tanh(par[1])
  total <- 0
  for (n in block_sizes) {
    C <- tridiag_corr(n, rho)
    total <- total + (-0.5) * determinant(C, logarithm = TRUE)$modulus[1]
  }
  total
}

## REML gradient
glm_weight_function <- function(mu, y, order_indices, family,
                                dispersion, observation_weights, ...) {
  rep(1, length(mu))
}
REML_grad <- function(par, model_fit, ...) {
  rho       <- tanh(par[1])
  drho_dpar <- 1 - rho^2   # sech^2(par[1])

  ## dV/drho is block-diagonal; each block is the tridiagonal matrix of 1s
  dV_drho_blocks <- lapply(block_sizes, function(n) {
    B <- matrix(0, n, n)
    if (n > 1) {
      for (i in 1:(n - 1)) {
        B[i, i + 1] <- 1
        B[i + 1, i] <- 1
      }
    }
    B
  })
  dV_drho <- collapse_block_diagonal(dV_drho_blocks)
  ## Chain rule: dV/dpar[1] = dV/drho * drho/dpar[1]
  dV_dpar1 <- drho_dpar * dV_drho
  grad1 <- reml_grad_from_dV(dV_dpar1, model_fit,
                             glm_weight_function, ...)
  c(grad1)
}

################################################################################

## ## Fitting ## ##
## Raw-scale analyses using identity link, assuming Gaussian response
model_fit <- lgspline(conc ~ spl(Time) + Time*Dose,
                data = df,
                K = unique(table(df$Subject)) - 1,
                correlation_id = df$Subject,
                VhalfInv_par_init = 0.1,
                VhalfInv_fxn = VhalfInv_fxn,
                Vhalf_fxn = Vhalf_fxn,
                VhalfInv_logdet = VhalfInv_logdet,
                REML_grad = REML_grad,
                spacetime = df$Time,
                qp_range_lower = 0,
                qp_observations =
                  list("Time:qp_negative_derivative" = which(df$Time > 5),
                       "Time:qp_positive_2ndderivative" = which(df$Time > 5),
                       "qp_range_lower" = which(df$Time>=0)),
                qp_negative_derivative = 'Time',
                qp_positive_2ndderivative = 'Time',
                include_cubic_terms = FALSE,
                include_quartic_terms = FALSE,
                include_constrain_second_deriv = FALSE,
                include_quadratic_interactions = TRUE,
                include_2way_interactions = FALSE)

## Closed-form equations for model fit
equation(model_fit)

################################################################################

## ## Non-Compartmental Analyses ## ##
#  Repeat for 1000 draws, posterior distributions of mean PK params for a subject
#  with a dose of 5 mg/L/kg
res <- sapply(1:1000, function(m){
  if(!(m %% 25)) cat('\n \t', m, ' / ', 1000)
  ## Generate posterior, including correlation, subject to constraints
  draw <- tryCatch(generate_posterior_correlation(
                             model_fit,
                             include_warnings = FALSE,
                             enforce_qp_constraints = FALSE,
                            ), error = function(e)
                            generate_posterior_correlation(
                              model_fit,
                              include_warnings = FALSE,
                              enforce_qp_constraints = FALSE
                            ))

  ## Find tmax and cmax
  find_tmax <- find_extremum(model_fit,
                             B_predict = draw$post_draw_coefficients,
                             sigmasq_predict = draw$post_draw_sigmasq,
                             initial = data.frame('Time'= 1.5,
                                                  'Dose'= 5))
  tmax <- find_tmax$t[,'Time']
  cmax <- find_tmax$y

  ## Find t_(1/2)
  find_thalf <- find_extremum(model_fit,
                              quick_heuristic = FALSE,
                              B_predict = draw$post_draw_coefficients,
                              minimize = TRUE,
                              initial = data.frame('Time' = 10,
                                                   'Dose' = 5),
                              custom_objective_function = function(mu,
                                                                   sigma,
                                                                   ybest,
                                                                   ...){
                                0.5*(mu - cmax/2)^2
                              },
                              custom_objective_derivative = function(mu,
                                                                     sigma,
                                                                     ybest,
                                                                     d_mu){
                                (mu - cmax/2) * d_mu
                              })
  thalf <- find_thalf$t[,'Time']

  ## Lambda-z based on tmax and thalf
  #  = (log(cmax/2) - log(cmax))/(tmax - thalf)
  lambdaz <- log(0.5)/(tmax - thalf)

  ## AUC-24
  auc24 <- integrate(model_fit,
                     B_predict = draw$post_draw_coefficients,
                     vars = 'Time',
                     lower = 0,
                     upper = 24)

  ## Predicted value at 24 hours
  yhat24 <- predict(model_fit,
                    B_predict = draw$post_draw_coefficients,
                    new_predictors = data.frame(
                      Time = 24,
                      Dose = 5
                    ))

  ## AUC-Inf
  aucInf <- auc24 + yhat24/lambdaz

  ## AUC-ratio
  aucratio <- auc24/aucInf

  ## Return all
  return(t(c(
    tmax,
    cmax,
    thalf,
    lambdaz,
    auc24,
    aucratio
  )))
})

## Filter for plausibility (AUC-24 > 0, AUC-Ratio < 1)
cat('\n\t ', ' Percent AUC-24<=0 or AUC-ratio >= 1: ', 100*mean(res[5,] <= 0 |
                                                                res[6,] >= 1),
    '% \n')
# > cat('\n\t ', ' Percent AUC-24<=0 or AUC-ratio >= 1: ', 100*mean(res[5,] <= 0 |
#                                                                     +                                                                 res[6,] >= 1),
#       +     '% \n')
#
# Percent AUC-24<=0 or AUC-ratio >= 1:  1.5 %
res <- res[,res[5,] > 0 & res[6,] < 1]

## Summarize results
for(i in 1:6){
  if(i == 1){
    title <- 'Time of Maximum Concentration (Hours)'
  }
  if(i == 2){
    title <- 'Maximum Concentration (mg/L)'
  }
  if(i == 3){
    title <- 'Time of Half-Maximum Concentration (mg/L)'
  }
  if(i == 4){
    title <- 'Terminal Rate Constant at 24 Hours'
  }
  if(i == 5){
    title <- 'AUC at 24 Hours'
  }
  if(i == 6){
    title <- 'AUC-24/AUC-Inf Ratio'
  }
  print(title)
  print(quantile(res[i,],
                 c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1),
                 na.rm = TRUE))
}
# [1] "Time of Maximum Concentration (Hours)"
# 0%      2.5%       25%       50%       75%     97.5%      100%
# 0.7799575 1.2685589 1.3374182 1.3822958 1.4581517 2.3121325 3.2589094
# [1] "Maximum Concentration (mg/L)"
# 0%      2.5%       25%       50%       75%     97.5%      100%
# 7.126560  7.841608  8.348087  8.609831  8.897443  9.491003 10.215822
# [1] "Time of Half-Maximum Concentration (mg/L)"
# 0%      2.5%       25%       50%       75%     97.5%      100%
# 7.120072  8.955549 10.371258 11.108551 11.758147 13.596257 24.650000
# [1] "Terminal Rate Constant at 24 Hours"
# 0%       2.5%        25%        50%        75%      97.5%       100%
# 0.02962187 0.05705144 0.06729804 0.07188852 0.07779173 0.09314024 0.12229288
# [1] "AUC at 24 Hours"
# 0%        2.5%         25%         50%         75%       97.5%        100%
# 0.8027022  41.4524128  77.8343274 100.4174735 122.1242632 163.2765239 200.6501400
# [1] "AUC-24/AUC-Inf Ratio"
# 0%       2.5%        25%        50%        75%      97.5%       100%
# 0.04700428 0.69161891 0.78611568 0.82648494 0.85932876 0.91797261 0.99144763

################################################################################

## ## More Inference
## Print out of individiual coefficients
Reduce("cbind", coef(model_fit))
cnf <- confint(model_fit)
print(cnf)

## Correlation estimate and crude 95% CI
print(tanh(model_fit$VhalfInv_params_est))
print(tanh(cnf[nrow(cnf),]))

## Wilk's test p-value for interaction term
#  Re-fit constraining interaction term to be 0
linear_constraint <- matrix(0, nrow = model_fit$P, ncol = model_fit$K + 1)
rownames(linear_constraint) <- names(unlist(model_fit$B))
for(k in 1:(model_fit$K + 1)){
  linear_constraint[paste0('partition', k, '.DosexTime^2'), k] <- 1
}
# Set their value to the 0 vector
null_value <- cbind(rep(0, ncol(linear_constraint)))
# Nested fit
nested_fit <-  lgspline(conc ~ spl(Time) + Time*Dose,
                        data = df,
                        K = unique(table(df$Subject)) - 1,
                        qp_observations = which(df$Time > 5),
                        qp_range_lower = 0,
                        qp_negative_derivative = 'Time',
                        qp_positive_2ndderivative = 'Time',
                        include_cubic_terms = FALSE,
                        include_quartic_terms = FALSE,
                        include_constrain_second_deriv = FALSE,
                        include_2way_interactions = FALSE,
                        include_quadratic_interactions = TRUE,
                        ## fixed args - same penalty, same covariance structure
                        #  same partitioning
                        opt = FALSE,
                        make_partition_list = model_fit$make_partition_list,
                        previously_tuned_penalties = model_fit$penalties,
                        VhalfInv = model_fit$VhalfInv,
                        Vhalf = Vhalf_fxn(model_fit$VhalfInv_params_estimates),
                        constraint_vectors = linear_constraint,
                        null_constraint = null_value
)

## Wilk's test
nll <- function(mod)sum(mod$VhalfInv %*% (mod$y - mod$ytilde)^2)*0.5/
                    model_fit$sigmasq_tilde
(stat <- -2*c(nll(model_fit) - nll(nested_fit)))
(degfed <- qr(nested_fit$A)$rank - qr(model_fit$A)$rank)
pchisq(
  q = stat,
  df = degfed,
  lower.tail = FALSE
)

# > nll <- function(mod)sum((mod$y - mod$ytilde)^2)*0.5/model_fit$sigmasq_tilde
# > (stat <- -2*c(nll(model_fit) - nll(nested_fit)))
# [1] 0.606065
# > (degfed <- qr(nested_fit$A)$rank - qr(model_fit$A)$rank)
# [1] 1
# > pchisq(
#   +   q = stat,
#   +   df = degfed,
#   +   lower.tail = FALSE
#   + )
# [1] 0.4362732
save.image(file = 'nca_results.RData')


################################################################################

## ## Two-Panel Figure: Fit + Posterior Densities ## ##
#  Left: fitted curve with posterior band over data
#  Right: 2x3 grid of posterior densities for the six PK summaries
load('nca_results.RData')
png("theoph_panel.png",
    width = 1800,
    height = 900,
    res = 120)

## Layout
layout(matrix(c(1, 1, 1, 2, 3, 4,
                1, 1, 1, 5, 6, 7),
              byrow = TRUE,
              nrow = 2))
## Left panel: fitted curve
par(mar = c(5, 5, 4, 2))
new_t <- seq(0, max(df$Time), length.out = 10000)
with(df, plot(Time, conc,
              main = 'Theophylline Concentration vs. Time',
              ylab = 'Concentration (mg / L)',
              xlab = 'Time (hours)',
              col = 'white',
              cex = 0,
              cex.axis = 1.4,
              cex.lab = 1.5,
              cex.main = 1.7,
              ylim = c(-0.1, 15)))
plot(model_fit,
     vars = 'Time',
     new_predictors = data.frame(
       'Time' = new_t,
       Dose = mean(df$Dose)
     ),
     add = TRUE,
     show_formulas = TRUE,
     include_all_terms_in_formulas = TRUE,
     digits = 4,
     cex = 1.4,
     custom_response_lab = 'Concentration',
     se.fit = TRUE,
     legend_args = list(cex = 1.2))
with(df, points(Time, conc, cex = 1.4))
## Right panels: posterior densities of the six PK summaries
pk_titles <- c(expression(t[max]~'(hours)'),
               expression(C[max]~'(mg/L)'),
               expression(t['1/2']~'(hours)'),
               expression(lambda[z]~'(h'^'-1'*')'),
               expression('AUC'['0-24']~'(mg'%.%'h/L)'),
               expression('AUC-24/AUC-'*infinity))
par(mar = c(4, 4, 3, 1))
for (i in 1:6) {
  dens <- density(res[i, ], na.rm = TRUE)
  q50 <- quantile(res[i, ], 0.50, na.rm = TRUE)
  q025 <- quantile(res[i, ], 0.025, na.rm = TRUE)
  q975 <- quantile(res[i, ], 0.975, na.rm = TRUE)

  plot(dens,
       main = pk_titles[i],
       xlab = '',
       ylab = 'Density',
       lwd = 2,
       col = 'black',
       cex.main = 1.3,
       cex.lab = 1.1,
       cex.axis = 1.0)

  ## Shade 95% credible interval
  ci_idx <- dens$x >= q025 & dens$x <= q975
  polygon(c(dens$x[ci_idx], rev(dens$x[ci_idx])),
          c(dens$y[ci_idx], rep(0, sum(ci_idx))),
          col = rgb(0.5, 0.5, 0.5, 0.4),
          border = NA)

  ## Median line
  abline(v = q50, lty = 2, lwd = 2, col = 'red')

  ## Redraw density on top
  lines(dens, lwd = 2)
}
dev.off()
