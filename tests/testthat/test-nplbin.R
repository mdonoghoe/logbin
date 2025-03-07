test_that("EM single observation works", {
  
  # Single observation intercept
  x_test <- matrix(1, ncol = 1)
  
  # y = 0 should throw an error
  expect_error(nplbin(y = 0, x = x_test, offset = 0, start = -1,
                      control = logbin.control(maxit = 1),
                      accelerate = "em"))
  
  # Also just test from 1:10
  N_test <- 10
  y_test <- seq(1, N_test - 1)
  
  for (i in seq_along(y_test)) {
    
    res <- nplbin(y = cbind(y_test[i], N_test - y_test[i]),
                  x = x_test, offset = 0, start = -1,
                  control = logbin.control(maxit = 1),
                  accelerate = "em")
    
    expect_equal(res$coefficients, log(y_test[i] / N_test))
    
  }
  
  # y = N should be at half the boundary
  boundary <- logbin.control()$bound.tol
  resN <- nplbin(y = cbind(N_test, 0),
                 x = x_test, offset = 0, start = -1,
                 control = logbin.control(maxit = 1),
                 accelerate = "em")
  expect_equal(resN$coefficients, -boundary / 2)

})

test_that("EM single step works", {
  
  x_test <- matrix(c(1,1,1,0,0.5,0.8), ncol = 2)
  y_test <- c(5, 6, 7)
  N_test <- rep(10, 3)
  beta_start <- c(-1, -4)
  
  # A simple implementation of the algorithm published by Marschner & Gillett 2012
  # NB we need the smallest non-zero x >= 1 for each column, so rescale x -> u first
  x_min <- apply(x_test, 2, function(t) min(t[t > 0]))
  u_test <- t(apply(x_test, 1, `/`, x_min))
  # (...u denotes a rescaled version of the parameter)
  betau_start <- x_min * beta_start
  lambdau_start <- exp(betau_start)
  p_start <- exp(u_test %*% betau_start)
  frac_start <- t(apply(p_start, 1, function(p) (lambdau_start - p) / (1 - p)))
  z_start <- apply(frac_start, 2, function(f) y_test + (N_test - y_test) * f)
  lambdau_update <- colSums(u_test * z_start) / (N_test %*% u_test)
  betau_update <- as.numeric(log(lambdau_update))
  # Put back on the original scale
  beta_update <- betau_update / x_min
  
  res <- nplbin(y = cbind(y_test, N_test - y_test),
                x = x_test, offset = 0, start = beta_start,
                control = logbin.control(maxit = 1),
                accelerate = "em")
  
  expect_equal(res$coefficients, beta_update)
  
})