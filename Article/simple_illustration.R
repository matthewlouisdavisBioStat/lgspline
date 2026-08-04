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

################################################################################

## ## GLM Extension ## ##
#  GLM extension (Section 3.2) applied to the same three-partition cubic
#  setup as above but with Bernoulli response.

## Penalty matrices
Lambda_s <- matrix(c(
  0, 0, 0,     0,
  0, 0, 0,     0,
  0, 0, 4 * 8, 6 * 8,
  0, 0, 6 * 8, 12 * 8
), byrow = TRUE, nrow = 4)
Lambda_r <- diag(c(1, 1, 0, 0)) 
Lambda <- 1e-4 * (Lambda_s + Lambda_r)
Lambda_Block <- Reduce("rbind", lapply(1:3, function(i){
  Reduce("cbind", lapply(1:3, function(j){
    Lambda * (j == i)
  }))
}))

## Tikhinov square root
LambdaHalf <- matsqrt(Lambda)

## Smoothness constraints, knots at -2 and 2
A <- t(matrix(c(
  1,  -2,   4,   -8,  -1,   2,  -4,    8,   0,   0,   0,    0,
  0,   1,  -4,   12,   0,  -1,   4,  -12,   0,   0,   0,    0,
  0,   0,   2,  -12,   0,   0,  -2,   12,   0,   0,   0,    0,
  0,   0,   0,    0,   1,   2,   4,    8,  -1,  -2,  -4,   -8,
  0,   0,   0,    0,   0,   1,   4,   12,   0,  -1,  -4,  -12,
  0,   0,   0,    0,   0,   0,   2,   12,   0,   0,  -2,  -12
), byrow = TRUE, nrow = 6))

## Create 3 partitions, knots at -2 and 2
t1 <- seq(-4, -2, length = 51)[1:50]
t2 <- seq(-2, 2, length = 51)[1:50]
t3 <- seq(2, 4, length = 50)
t <- c(t1, t2, t3)

## Create response
logit_probs <- cos(t) +
  -1 +
  100 * t +
  0.025 * t^2 +
  (1/60) * t^3 +
  (1/100) * t^5 -
  0.0025 * t^6 +
  0.01 * (exp(t) - 1) +
  log(abs(t) + 1) +
  5 * sin(t) +
  1 / (1 + exp(-t)) -
  0.05 * t^4

## Standardize and sigmoid transform
probs <- 0.5 * (plogis(std(logit_probs)) + 0.5)

## Generate outcome
y <- sapply(probs, function(p)rbinom(1, 1, p))

## Partition
y1 <- y[1:50]
y2 <- y[1:50 + 50]
y3 <- y[1:50 + 100]

## Design matrices per partition
X1 <- cbind(1, t1, t1^2, t1^3)
X2 <- cbind(1, t2, t2^2, t2^3)
X3 <- cbind(1, t3, t3^2, t3^3)

## Block-diagonal full design
X <- rbind(cbind(X1,   0*X1, 0*X1),
           cbind(0*X2, X2,   0*X2),
           cbind(0*X3, 0*X3, X3))

## Tikhinov-augmented per-partition design for the GLM warm start
X1_fit <- rbind(X1, LambdaHalf)
X2_fit <- rbind(X2, LambdaHalf)
X3_fit <- rbind(X3, LambdaHalf)
y1_fit <- c(y1, rep(quasibinomial()$linkinv(0), nrow(LambdaHalf)))
y2_fit <- c(y2, rep(quasibinomial()$linkinv(0), nrow(LambdaHalf)))
y3_fit <- c(y3, rep(quasibinomial()$linkinv(0), nrow(LambdaHalf)))

## Unconstrained penalized MAP estimators per partition
fit1 <- glm.fit(X1_fit, y1_fit, family = quasibinomial())
fit2 <- glm.fit(X2_fit, y2_fit, family = quasibinomial())
fit3 <- glm.fit(X3_fit, y3_fit, family = quasibinomial())
bhat1 <- fit1$coef
bhat2 <- fit2$coef
bhat3 <- fit3$coef

## Clean up with weighted Newton-Raphson steps
yhat1 <- plogis(X1 %*% bhat1)
yhat2 <- plogis(X2 %*% bhat2)
yhat3 <- plogis(X3 %*% bhat3)
for(i in 1:10){
  W1 <- diag(c(yhat1 * (1 - yhat1)))
  W2 <- diag(c(yhat2 * (1 - yhat2)))
  W3 <- diag(c(yhat3 * (1 - yhat3)))
  bhat1 <- 0.9 * bhat1 +
    0.1 * solve(t(X1) %*% W1 %*% X1 + Lambda) %*% t(X1) %*% (y1 - yhat1)
  bhat2 <- 0.9 * bhat2 +
    0.1 * solve(t(X2) %*% W2 %*% X2 + Lambda) %*% t(X2) %*% (y2 - yhat2)
  bhat3 <- 0.9 * bhat3 +
    0.1 * solve(t(X3) %*% W3 %*% X3 + Lambda) %*% t(X3) %*% (y3 - yhat3)
  yhat1 <- plogis(X1 %*% bhat1)
  yhat2 <- plogis(X2 %*% bhat2)
  yhat3 <- plogis(X3 %*% bhat3)
}
bhat <- c(bhat1, bhat2, bhat3)

## Apply smoothness projection via transformed-OLS, iterating G(beta)
btilde <- bhat
for(i in 1:5){
  mu <- plogis(X %*% btilde)
  W  <- diag(c(mu * (1 - mu)))
  G  <- solve(t(X) %*% W %*% X + Lambda_Block)
  Ghalf    <- matsqrt(G)
  GhalfInv <- matinvsqrt(G)
  ystar <- GhalfInv %*% btilde
  Xstar <- Ghalf %*% A
  gamma <- solve(t(Xstar) %*% Xstar, t(Xstar) %*% ystar)
  btilde <- as.numeric(Ghalf %*% (ystar - Xstar %*% gamma))
}

## Fitted values
ytilde1 <- plogis(X1 %*% btilde[1:4])
ytilde2 <- plogis(X2 %*% btilde[5:8])
ytilde3 <- plogis(X3 %*% btilde[9:12])

## Equation labels
eq1 <- sprintf("logit(p) = %.2f + %.2ft + %.2ft^2 + %.2ft^3",
               btilde[1], btilde[2], btilde[3], btilde[4])
eq2 <- sprintf("logit(p) = %.2f + %.2ft + %.2ft^2 + %.2ft^3",
               btilde[5], btilde[6], btilde[7], btilde[8])
eq3 <- sprintf("logit(p) = %.2f + %.2ft + %.2ft^2 + %.2ft^3",
               btilde[9], btilde[10], btilde[11], btilde[12])

## Plot fitted probabilities
plot(t, y,
     main = "GLM Fit with Smoothness Constraints",
     cex = 2.5, cex.lab = 2, cex.main = 2,
     ylim = c(-0.195, 1), xlab = "t")
points(t1, ytilde1, col = 'darkblue',   type = 'l', lwd = 8)
points(t2, ytilde2, col = 'blue',       type = 'l', lwd = 8)
points(t3, ytilde3, col = 'darkviolet', type = 'l', lwd = 8)
legend("bottomright",
       legend = c(eq1, eq2, eq3),
       col = c("darkblue", "blue", "darkviolet"),
       lwd = 4, cex = 1.2, bg = "white", inset = c(0.02))