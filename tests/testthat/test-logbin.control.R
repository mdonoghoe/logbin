test_that("logbin.control requires positive boundary tolerance", {
  
  expect_error(logbin.control(bound.tol = "a"))
  
  expect_error(logbin.control(bound.tol = -100))
  
  expect_error(logbin.control(bound.tol = 0))
  
})

test_that("logbin.control requires positive epsilon", {
  
  expect_error(logbin.control(epsilon = "a"))
  
  expect_error(logbin.control(epsilon = -100))
  
  expect_error(logbin.control(epsilon = 0))
  
})

test_that("logbin.control expects epsilon < bound.tol", {
  
  expect_warning(logbin.control(epsilon = 1e-8, bound.tol = 1e-8))
  
  expect_warning(logbin.control(epsilon = 1e-6, bound.tol = 1e-8))
  
})

test_that("logbin.control requires positive maxit", {
  
  expect_error(logbin.control(maxit = "a"))
  
  expect_error(logbin.control(maxit = -100))
  
  expect_error(logbin.control(maxit = 0))
  
})