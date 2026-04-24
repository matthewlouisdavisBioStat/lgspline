test_that("basic plotting functions work", {
  set.seed(1234)
  t <- seq(-9, 9, length.out = 100)
  y <- sin(t) + rnorm(100, 0, 0.1)
  fit <- lgspline(cbind(t), y, K = 5)
  # Test 1D plot
  expect_error(plot(fit), NA)
  expect_error(plot(fit, show_formulas = TRUE), NA)
  # Test 2D plot, include quartic terms
  data(volcano)
  volcano_long <- Reduce('rbind', lapply(1:nrow(volcano), function(i){
    t(sapply(1:ncol(volcano), function(j){
      c(i, j, volcano[i,j])
    }))
  }))
  fit2d <- lgspline(volcano_long[,1:2],
                    include_quartic_terms = TRUE,
                    volcano_long[,3], K = 5)
  plot(fit2d,
       show_formulas = TRUE,
       text_size_formula  = 7,
       digits = 2,
       custom_predictor_lab1 = 'My Predictor 1',
       custom_predictor_lab2 = 'My Predictor 2',
       custom_response_lab = 'My Response')
})
