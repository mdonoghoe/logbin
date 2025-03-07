test_that("predict.logbin gives same results as glm for a simple model", {
  
  require(glm2)
  data(heart)
  fit.glm <- glm(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup),
                 family = binomial(log), data = heart,
                 start = rep(-0.5, 3))
  
  fit.logbin <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup), 
                       data = heart, start = rep(-0.5, 3), method = "glm")
  
  expect_equal(predict(fit.glm, newdata = NULL, type = "link", 
                       se.fit = FALSE, na.action = na.pass),
               predict(fit.logbin, newdata = NULL, type = "link",
                       na.action = na.pass))
  
  expect_equal(predict(fit.glm, newdata = NULL, type = "response",
                       se.fit = FALSE, na.action = na.pass),
               predict(fit.logbin, newdata = NULL, type = "response",
                       na.action = na.pass))
  
})

test_that("predict.logbin with type = terms gives same results as glm for a simple model", {
  
  require(glm2)
  data(heart)
  fit.glm <- glm(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup),
                 family = binomial(log), data = heart,
                 start = rep(-0.5, 3))
  
  fit.logbin <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup), 
                       data = heart, start = rep(-0.5, 3), method = "glm")
  
  pred.glm <- predict(fit.glm, newdata = NULL, type = "terms",
                      se.fit = FALSE, na.action = na.pass)
  
  pred.logbin <- predict(fit.logbin, newdata = NULL, type = "terms",
                         na.action = na.pass)
  
  expect_equal(pred.glm, pred.logbin)
  expect_equal(attr(pred.glm, "constant"), attr(pred.logbin, "constant"))
  expect_identical(dim(pred.glm), dim(pred.logbin))
  expect_identical(colnames(pred.glm), colnames(pred.logbin))
  
})

test_that("predict.logbin gives same results as glm with newdata", {
  
  require(glm2)
  data(heart)
  
  # Fit Region as a continuous variable so we can use other values
  fit.glm <- glm(cbind(Deaths, Patients-Deaths) ~ Region,
                 family = binomial(log), data = heart,
                 start = rep(-0.5, 2))
  
  fit.logbin <- logbin(cbind(Deaths, Patients-Deaths) ~ Region, 
                       data = heart, start = rep(-0.5, 2), method = "glm")
  
  newdat <- data.frame(Region = c(seq(1, 3, 0.1), NA))
  
  expect_equal(predict(fit.glm, newdata = newdat, type = "link",
                       se.fit = FALSE, na.action = na.pass),
               predict(fit.logbin, newdata = newdat, type = "link",
                       se.fit = FALSE, na.action = na.pass))
  
  expect_equal(predict(fit.glm, newdata = newdat, type = "link",
                       se.fit = FALSE, na.action = na.omit),
               predict(fit.logbin, newdata = newdat, type = "link",
                       se.fit = FALSE, na.action = na.omit))
  
  expect_equal(predict(fit.glm, newdata = newdat, type = "response",
                       se.fit = FALSE, na.action = na.omit),
               predict(fit.logbin, newdata = newdat, type = "response",
                       se.fit = FALSE, na.action = na.omit))
  
  # Should error due to the NA
  expect_error(predict(fit.logbin, newdata = newdat, type = "response",
                       na.action = na.fail))
  
})

test_that("predict.logbin checks xminmax", {
  
  require(glm2)
  data(heart)
  
  fit.logbin <- logbin(cbind(Deaths, Patients-Deaths) ~ Region, 
                       data = heart, start = rep(-0.5, 2), method = "em")
  
  newdat <- data.frame(Region = c(seq(1, 3.5, 0.1), NA))
  newdat_ok <- subset(newdat, Region <= 3)
  
  # Should error trying to predict for Region > 3
  expect_error(predict(fit.logbin, newdata = newdat, type = "response",
                       checkminmax = TRUE))
  
  # No error if we ignore the check...
  expect_vector(predict(fit.logbin, newdata = newdat, type = "response",
                        checkminmax = FALSE), 
                ptype = double(), size = nrow(newdat))
  # ...or if we remove the values > 3
  expect_vector(predict(fit.logbin, newdata = newdat_ok, type = "response",
                        checkminmax = TRUE),
                ptype = double(), size = nrow(newdat_ok))
  
  
})