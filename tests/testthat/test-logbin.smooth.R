test_that("logbin.smooth with B and fixed knots works", {
  
  data(mtcars)
  
  suppressWarnings(mdl <- logbin.smooth(am ~ B(mpg, knots = c(15, 20, 25)), 
                                        data = mtcars, method = "em", maxit = 1e10))
  
  expect_s3_class(mdl, "logbin.smooth")
  
  expect_true(all(predict(mdl) < 0))
  
})

test_that("logbin.smooth with B and knot range works", {
  
  data(mtcars)
  
  suppressWarnings(mdl <- logbin.smooth(am ~ B(mpg, knot.range = 0:5), 
                                        data = mtcars, method = "em", maxit = 1e10))
  
  expect_s3_class(mdl, "logbin.smooth")
  
  expect_true(all(predict(mdl) < 0))

})

test_that("logbin.smooth with B and monotonic restriction works", {

  data(mtcars)
  
  suppressWarnings(mdl <- logbin.smooth(am ~ B(mpg, knots = c(15, 20, 25)), mono = 1,
                                        data = mtcars, method = "em", maxit = 1e10))
  
  expect_s3_class(mdl, "logbin.smooth")
  
  expect_true(all(predict(mdl) < 0))
  
  # Check that predictions are monotonic
  expect_true(all(diff(predict(mdl)[order(mtcars$mpg)]) >= 0))

})

test_that("logbin.smooth with Iso works", {
  
  data(mtcars)
  
  suppressWarnings(mdl <- logbin.smooth(am ~ Iso(mpg), mono = 1,
                                        data = mtcars, method = "em", maxit = 1e10))
  
  expect_s3_class(mdl, "logbin.smooth")
  
  expect_true(all(predict(mdl) < 0))
  
  # Check that predictions are monotonic
  expect_true(all(diff(predict(mdl)[order(mtcars$mpg)]) >= 0))

})