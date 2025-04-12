test_that("basic correlation structure runs without error", {
  set.seed(1234)
  x <- seq(-9, 9, length.out = 1000)
  y <- sin(x) + rnorm(1000, 0, 0.1)
  matsqrt <- function(mat){
    eig <- eigen(mat)
    eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  }
  model_fit <- lgspline(x, y, K = 1, standardize_response = FALSE,
                  VhalfInv = diag(length(y)))
  expect_error(model_fit$model_fit$VhalfInv_params_estimates, NA)
})

