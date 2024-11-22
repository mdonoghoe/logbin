test_that("EM single step works", {
  # Do this randomly so we don't fluke a pass?
  N_test <- 2 + rbinom(1, size = 48, prob = 0.8)
  # Ensure we don't get zero (-Inf)
  y_test <- 1 + rbinom(1, size = N_test - 1, prob = 0.3)
  
  res <- nplbin(y = cbind(y_test, N_test - y_test), 
                x = matrix(1, ncol = 1),
                offset = 0, start = -1,
                control = logbin.control(maxit = 1),
                accelerate = "em")

  expect_equal(res$coefficients, log(y_test / N_test))
})