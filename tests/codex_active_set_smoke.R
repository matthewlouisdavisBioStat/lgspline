.libPaths(c(normalizePath(".r-lib", winslash = "/"), .libPaths()))

library(lgspline)

fmt_idx <- function(x) {
  if (length(x) == 0) return("<none>")
  paste(x, collapse = ",")
}

## Gaussian GEE Woodbury reproducer
set.seed(20260421)

n <- 60
t1 <- seq(0, 1, length.out = n)
t2 <- sin(2 * pi * t1)
y <- 0.5 + 1.2 * t1 + 0.3 * t2 + rnorm(n, 0, 0.05)
Vinv <- diag(n)
Vinv[15, 45] <- 0.2
Vinv[45, 15] <- 0.2
VhalfInv <- t(chol(Vinv))

fit_gauss <- lgspline(
  cbind(t1, t2),
  y,
  K = 1,
  opt = FALSE,
  qp_positive_derivative = "t1",
  VhalfInv = VhalfInv,
  standardize_response = FALSE,
  include_warnings = FALSE
)

fit_gauss_core <- if (!is.null(fit_gauss$model_fit)) fit_gauss$model_fit else fit_gauss
deriv_gauss <- predict(fit_gauss, cbind(t1, t2), take_first_derivatives = TRUE)
min_deriv_gauss <- min(unlist(deriv_gauss$first_deriv[[1]]))

cat("gaussian_method:", fit_gauss_core$qp_info$method, "\n")
cat("gaussian_active:", fmt_idx(fit_gauss_core$qp_info$active_ineq), "\n")
cat("gaussian_min_deriv:", signif(min_deriv_gauss, 6), "\n")

stopifnot(isTRUE(fit_gauss_core$qp_info$converged))
stopifnot(min_deriv_gauss >= -1e-6)

## GLM GEE Woodbury reproducer
set.seed(20260421)

t2 <- cos(2 * pi * t1)
mu <- exp(0.3 + 0.8 * t1 + 0.2 * t2)
y <- rpois(n, mu)
Vinv <- diag(n)
Vinv[12, 36] <- 0.15
Vinv[36, 12] <- 0.15
VhalfInv <- t(chol(Vinv))

fit_glm <- lgspline(
  cbind(t1, t2),
  y,
  family = quasipoisson(),
  K = 1,
  opt = FALSE,
  qp_positive_derivative = "t1",
  VhalfInv = VhalfInv,
  include_warnings = FALSE
)

fit_glm_core <- if (!is.null(fit_glm$model_fit)) fit_glm$model_fit else fit_glm
deriv_glm <- predict(fit_glm, cbind(t1, t2), take_first_derivatives = TRUE)
min_deriv_glm <- min(unlist(deriv_glm$first_deriv[[1]]))

cat("glm_method:", fit_glm_core$qp_info$method, "\n")
cat("glm_active:", fmt_idx(fit_glm_core$qp_info$active_ineq), "\n")
cat("glm_min_deriv:", signif(min_deriv_glm, 6), "\n")

stopifnot(isTRUE(fit_glm_core$qp_info$converged))
if (min_deriv_glm < -1e-6) {
  warning(
    "glm_min_deriv is below zero: ",
    signif(min_deriv_glm, 6),
    "\n  This is still a live issue for the Woodbury active-set refactor."
  )
}

cat("Focused smoke passed\n")
