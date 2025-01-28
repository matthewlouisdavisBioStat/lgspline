test_that("lgspline works with parallel processing", {
  skip_on_cran()  # Skip on CRAN out of common courtesy

  # Setup test data
  set.seed(1234)
  x <- seq(-9, 9, length.out = 1000)
  y <- sin(x) + rnorm(1000, 0, 0.1)
  dat <- cbind(x, y)

  # Test parallel vs non-parallel results match
  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl))
  ## Ensure cluster is stopped even if test fails

  set.seed(1234)
  fit_parallel <- lgspline(cbind(dat[,'x']),
                          dat[,'y'],
                          cl = cl,
                          K = 2)

  set.seed(1234)
  fit_serial <- lgspline(cbind(dat[,'x']),
                        dat[,'y'],
                        K = 2)

  # Compare results
  expect_equal(fit_parallel$ytilde,
              fit_serial$ytilde,
              tolerance = 1e-5)

  # Test predictions match
  newx <- matrix(seq(-9, 9, length.out = 100))
  pred_parallel <- predict(fit_parallel, newx)
  pred_serial <- predict(fit_serial, newx)

  expect_equal(pred_parallel, pred_serial, tolerance = 1e-5)
})

test_that("lgspline parallel processing handles errors and options", {
  skip_on_cran() # Skip on CRAN out of common courtesy

  set.seed(1234)
  x1 <- seq(-9, 9, length.out = 1000)
  x2 <- seq(-9, 9, length.out = 1000)
  y <- sin(x1) + cos(x2) + rnorm(1000, 0, 0.1)

  # Test with different parallel options
  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl))

  # Test aga penalties, turn off eigen penalties
  expect_no_error(
    lgspline(cbind(x1,x2),
             y,
             cl = cl,
             K = 1,
             parallel_eigen = FALSE,
             parallel_aga = TRUE,
             include_warnings = FALSE)
  )

  # Test parallel penalty/constraint/neighbor computation
  expect_no_error(
    lgspline(predictors = cbind(x1,x2),
             response = y,
             cl = cl,
             K = 2,
             parallel_find_neighbors = TRUE,
             parallel_trace = TRUE,
             parallel_matmult = TRUE,
             parallel_make_constraint = TRUE,
             parallel_penalty = TRUE,
             include_warnings = FALSE)
  )
})
