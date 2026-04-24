test_that("basic lgspline fits 1D data correctly using base R BFGS for tuning", {
  t <- seq(-9, 9, length.out = 100)
  y <- sin(t) + rnorm(100, 0, 0.1)
  fit <- lgspline(cbind(t), y, K = 2, use_custom_bfgs = FALSE)
  expect_s3_class(fit, "lgspline")
  expect_length(fit$B, fit$K + 1)
  expect_true(all(!is.na(fit$ytilde)))
})
test_that("Basic lgspline fits 2D data correctly", {
  data(volcano)
  volcano_long <- Reduce('rbind', lapply(1:nrow(volcano), function(i){
    t(sapply(1:ncol(volcano), function(j){
      c(i, j, volcano[i,j])
    }))
  }))
  fit <- lgspline(volcano_long[,1:2], volcano_long[,3], K = 1, opt = FALSE,
                  unique_penalty_per_partition = FALSE)
  expect_s3_class(fit, "lgspline")
  expect_true(fit$q == 2)
  expect_true(all(!is.na(fit$ytilde)))
})
test_that("Basic lgspline handles no knots", {
  # No knots
  t <- 1:10
  y <- t + rnorm(10)
  expect_error(lgspline(cbind(t), y, K = 0), NA)
})
