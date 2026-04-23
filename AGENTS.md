# Codex Notes for lgspline

## Start here

1. Run `Rscript scripts/codex_setup.R`
2. Reinstall the package using the PowerShell command printed by setup
3. Read:
   - `R/solver_utils.R`
   - `R/get_B.R`
   - `R/blockfit_solve.R`
   - `tests/testthat/test_correlation_structure.R`
   - `tests/codex_active_set_smoke.R`
4. Run `Rscript scripts/codex_verify.R`

## Current goal

Use one generic active-set wrapper around equality-only solves.
Do not gate active-set on partition-local vs dense structure.
Dense QP/SQP should stay fallback-only.

## Current status

- Shared active-set helpers now live in `R/solver_utils.R`
- `get_B()` and `blockfit_solve()` now try active-set first
- GLM Woodbury weight blowups were stabilized enough to avoid the old crash
- The two focused reproductions still fall back:
  - Gaussian GEE Woodbury -> `dense_qp_gee_gaussian`
  - GLM GEE Woodbury -> `dense_qp_gee_glm`
- Gaussian fallback gives a sparse active set and near-nonnegative derivative
- GLM fallback still has negative derivative mass in the focused smoke
- The missing piece is making the Woodbury equality re-solve / KKT check accept that active set

## Useful facts

- Prefer reinstalling the package into `.r-lib` from PowerShell
- `pkgload::load_all()` was unreliable here
- `source()` alone is not enough because the package uses compiled code
- `tests/codex_active_set_smoke.R` is the quickest reproducer

## Comment style

- Keep comments short
- Prefer `## Section` and, if needed, `#   detail`, see example below
- align notation and "flow" with the overleaf document when presenting comments
- keep code efficient, try to re-use and copy my style as much as possible while
  keeping code correct, compatible, and goal-achieving without sidestepping or
  taking shortcuts.
- Raw predictors should be named `t`, `t1`, `t2`, etc., never `x1` or `x2`
- Reserve `x` names for expansions or design matrices

set.seed(1234)
library(lgspline)

## Penalty matrices, using the cubic smoothing spline penalty
Lambda_s <- matrix(c(
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 4 * (4 - -4), 6 * (4 - -4),
  0, 0, 6 * (4 - -4), 12 * (4 - -4)
), byrow = TRUE, nrow = 4)
Lambda_r <- diag(4) / 25 # flat ridge penalty
Lambda_r[-c(1:2),-c(1:2)] <- 0 # quadratic and cubic terms have no ridge
Lambda <- 0.01 * (Lambda_r + Lambda_s) # lambda_w = 0.01
Lambda <- Reduce("rbind", lapply(1:3, function(i){
  Reduce("cbind", lapply(1:3, function(j){
    Lambda * (j == i)
  }))
}))

## Create 3 partitions, knots at -2 and 2
t1 <- seq(-4, -2, length = 51)[1:50]
t2 <- seq(-2, 2, length = 51)[1:50]
t3 <- seq(2, 4, length = 50)
t <- c(t1, t2, t3)

## Create response, and partition
y <- -1 +
  0.2 * t +
  0.025 * t^2 +
  (1/60) * t^3 +
  (1/100) * t^5 -
  0.0025 * t^6 +
  0.01 * (exp(t) - 1) +
  log(abs(t) + 1) +
  5 * sin(t) +
  1/(1 + exp(-t)) -
  0.05 * t^4 +
  rnorm(length(t), 0, 10)
y1 <- y[1:50]
y2 <- y[1:50 + 50]
y3 <- y[1:50 + 100]

## Create design matrices
X1 <- cbind(1, t1, t1^2, t1^3)
X2 <- cbind(1, t2, t2^2, t2^3)
X3 <- cbind(1, t3, t3^2, t3^3)

## Get coefficients (Step 1: unconstrained estimates)
bhat1 <- solve(t(X1) %*% X1 + Lambda[1:4, 1:4]) %*% t(X1) %*% y1
bhat2 <- solve(t(X2) %*% X2 + Lambda[5:8, 5:8]) %*% t(X2) %*% y2
bhat3 <- solve(t(X3) %*% X3 + Lambda[9:12, 9:12]) %*% t(X3) %*% y3
bhat <- c(bhat1, bhat2, bhat3)

## Fitted values (unconstrained)
yhat1 <- X1 %*% bhat1
yhat2 <- X2 %*% bhat2
yhat3 <- X3 %*% bhat3

## Computing G partition by partition (stored as list)
G_list <- list()
G_list[[1]] <- solve(t(X1) %*% X1 + Lambda[1:4, 1:4])
G_list[[2]] <- solve(t(X2) %*% X2 + Lambda[5:8, 5:8])
G_list[[3]] <- solve(t(X3) %*% X3 + Lambda[9:12, 9:12])

## Code smoothing constraints
A <- t(matrix(c(
  1,  -2,   4,   -8,  -1,   2,  -4,    8,   0,   0,   0,    0, # fitted at -2
  0,   1,  -4,   12,   0,  -1,   4,  -12,   0,   0,   0,    0, # first deriv at -2,
  0,   0,   2,  -12,   0,   0,  -2,   12,   0,   0,   0,    0, # second deriv at -2,
  0,   0,   0,    0,   1,   2,   4,    8,  -1,  -2,  -4,   -8, # fitted at 2
  0,   0,   0,    0,   0,   1,   4,   12,   0,  -1,  -4,  -12, # first deriv at 2,
  0,   0,   0,    0,   0,   0,   2,   12,   0,   0,  -2,  -12  # second deriv at 2
), byrow = TRUE, nrow = 6)
)

## Compute G^(1/2) and G^(-1/2) partition by partition via eigendecomposition
G_sqrt_list <- list()
G_inv_sqrt_list <- list()

for(k in 1:3) {
  G_sqrt_list[[k]] <- matsqrt(G_list[[k]])
  G_inv_sqrt_list[[k]] <- matinvsqrt(G_list[[k]])
}

## Construct block-diagonal G^(1/2) and G^(-1/2) from partition-specific matrices
G_sqrt <- matrix(0, nrow = 12, ncol = 12)
G_inv_sqrt <- matrix(0, nrow = 12, ncol = 12)

G_sqrt[1:4, 1:4] <- G_sqrt_list[[1]]
G_sqrt[5:8, 5:8] <- G_sqrt_list[[2]]
G_sqrt[9:12, 9:12] <- G_sqrt_list[[3]]

G_inv_sqrt[1:4, 1:4] <- G_inv_sqrt_list[[1]]
G_inv_sqrt[5:8, 5:8] <- G_inv_sqrt_list[[2]]
G_inv_sqrt[9:12, 9:12] <- G_inv_sqrt_list[[3]]

## Transformed OLS
#  Step 2: Transform y* = G^(-1/2) * bhat and X* = G^(1/2) * A
y_star <- G_inv_sqrt %*% bhat
X_star <- G_sqrt %*% A
#  Step 3: Fit OLS model E[y*] = X* gamma*
beta_star <- solve(t(X_star) %*% X_star, t(X_star) %*% y_star)
#  Step 4: Obtain residuals r* and compute constrained estimate
r_star <- y_star - X_star %*% beta_star
btilde <- G_sqrt %*% r_star

## Constrained fitted values
ytilde1 <- X1 %*% btilde[1:4]
ytilde2 <- X2 %*% btilde[5:8]
ytilde3 <- X3 %*% btilde[9:12]

## Create legend labels with equations
# unconstrained, by partition
eq1 <- sprintf("f(t) = %.2f + %.2ft + %.2ft^2 + %.2ft^3",
               bhat1[1], bhat1[2], bhat1[3], bhat1[4])
eq2 <- sprintf("f(t) = %.2f + %.2ft + %.2ft^2 + %.2ft^3",
               bhat2[1], bhat2[2], bhat2[3], bhat2[4])
eq3 <- sprintf("f(t) = %.2f + %.2ft + %.2ft^2 + %.2ft^3",
               bhat3[1], bhat3[2], bhat3[3], bhat3[4])

# constrained, by partition
eq1_c <- sprintf("f(t) = %.2f + %.2ft + %.2ft^2 + %.2ft^3",
                 btilde[1], btilde[2], btilde[3], btilde[4])
eq2_c <- sprintf("f(t) = %.2f + %.2ft + %.2ft^2 + %.2ft^3",
                 btilde[5], btilde[6], btilde[7], btilde[8])
eq3_c <- sprintf("f(t) = %.2f + %.2ft + %.2ft^2 + %.2ft^3",
                 btilde[9], btilde[10], btilde[11], btilde[12])

png("simple_illustration.png",
    width = 1800,
    height = 900,
    res = 120)

## Plot results with legends
layout(matrix(c(1, 1, 1, 2, 2, 2,
                1, 1, 1, 2, 2, 2,
                1, 1, 1, 2, 2, 2),
              nrow = 3,
              byrow = TRUE))

## First plot (Unconstrained)
plot(t, y, main = "Unconstrained Fit", cex = 2.5, cex.lab = 2, cex.main = 2, xlab = "t")
points(t1, yhat1, col = 'darkblue', type = 'l', lwd = 8)
points(t2, yhat2, col = 'blue', type = 'l', lwd = 8)
points(t3, yhat3, col = 'darkviolet', type = 'l', lwd = 8)
legend("bottomright",
       legend = c(eq1, eq2, eq3),
       col = c("darkblue", "blue", "darkviolet"),
       lwd = 4,
       cex = 2,
       bg = "white",
       inset = 0.02)

## Second plot (Constrained)
plot(t, y, main = "Fit Under Smoothing Constraints",
     cex = 2.5, cex.lab = 2, cex.main = 2, xlab = "t")
points(t1, ytilde1, col = 'darkblue', type = 'l', lwd = 8)
points(t2, ytilde2, col = 'blue', type = 'l', lwd = 8)
points(t3, ytilde3, col = 'darkviolet', type = 'l', lwd = 8)
legend("bottomright",
       legend = c(eq1_c, eq2_c, eq3_c),
       col = c("darkblue", "blue", "darkviolet"),
       lwd = 4,
       cex = 2,
       bg = "white",
       inset = 0.02)

dev.off()
