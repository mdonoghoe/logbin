test_that("logbin gives same results as glm for a simple model", {
  
  require(glm2)
  data(heart)
  fit.glm <- glm(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup),
                 family = binomial(log), data = heart,
                 start = rep(-0.5, 3))
  
  fit.logbin <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup), 
                       data = heart, start = rep(-0.5, 3), method = "glm")
  
  common_parts <- c("coefficients", "residuals", "fitted.values")
  
  expect_equal(fit.glm[common_parts], fit.logbin[common_parts])
  
})

test_that("logbin finds a higher likelihood than glm in heart data", {
  
  require(glm2)
  data(heart)
  
  start.p <- sum(heart$Deaths) / sum(heart$Patients)
  
  # glm will throw warnings due to step size
  suppressWarnings(
    fit.glm <- glm(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup) +
                     factor(Severity) + factor(Delay) + factor(Region),
                   family = binomial(log), data = heart,
                   start = c(log(start.p), rep(c(0.2, 0.4), 4)),
                   maxit = 100)
  )
  
  fit.logbin <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup) +
                         factor(Severity) + factor(Delay) + factor(Region),
                       data = heart,
                       start = c(log(start.p), rep(c(0.2, 0.4), 4)))
  
  expect_gte(fit.logbin$loglik, logLik(fit.glm))
  
})

test_that("accelerate options produce the same results", {
  
  require(glm2)
  data(heart)
  
  start.p <- sum(heart$Deaths) / sum(heart$Patients)
  
  # Fit the model with the default (EM)...
  fit.logbin <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup) +
                         factor(Severity) + factor(Delay) + factor(Region),
                       data = heart,
                       start = c(log(start.p), rep(c(0.2, 0.4), 4)))
  
  # ... and the accelerated EM algorithms
  fit.logbin.sqem <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup) +
                              factor(Severity) + factor(Delay) + factor(Region),
                            data = heart, accelerate = "squarem",
                            start = c(log(start.p), rep(c(0.2, 0.4), 4)))
  
  fit.logbin.pem <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup) +
                             factor(Severity) + factor(Delay) + factor(Region),
                           data = heart, accelerate = "pem",
                           start = c(log(start.p), rep(c(0.2, 0.4), 4)))
  
  fit.logbin.qn <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup) +
                            factor(Severity) + factor(Delay) + factor(Region),
                          data = heart, accelerate = "qn",
                          start = c(log(start.p), rep(c(0.2, 0.4), 4)))
  
  # Log-likelihood should be almost identical...
  expect_equal(fit.logbin$loglik, fit.logbin.sqem$loglik)
  expect_equal(fit.logbin$loglik, fit.logbin.pem$loglik)
  expect_equal(fit.logbin$loglik, fit.logbin.qn$loglik)
  
  # Allow a bit more tolerance in the coefficients
  expect_equal(fit.logbin$coefficients, fit.logbin.sqem$coefficients,
               tolerance = 1e-4)
  expect_equal(fit.logbin$coefficients, fit.logbin.pem$coefficients,
               tolerance = 1e-4)
  expect_equal(fit.logbin$coefficients, fit.logbin.qn$coefficients,
               tolerance = 1e-4)
  
})

test_that("method options produce the same results", {
  
  require(glm2)
  data(heart)
  
  # Simple model so that we can compare with glm/glm2...
  # ...which don't work with the more full model
  
  fit.logbin <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup), 
                       data = heart, start = rep(-0.5, 3))
  
  fit.logbin.em <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup), 
                          data = heart, start = rep(-0.5, 3), method = "em")
  
  fit.logbin.glm <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup), 
                           data = heart, start = rep(-0.5, 3), method = "glm")
  
  fit.logbin.glm2 <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup), 
                            data = heart, start = rep(-0.5, 3), method = "glm2")
  
  fit.logbin.ab <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup), 
                          data = heart, start = rep(-0.5, 3), method = "ab")
  
  # Log-likelihood should be almost identical...
  expect_equal(fit.logbin$loglik, fit.logbin.em$loglik)
  expect_equal(fit.logbin$loglik, fit.logbin.glm$loglik)
  expect_equal(fit.logbin$loglik, fit.logbin.glm2$loglik)
  expect_equal(fit.logbin$loglik, fit.logbin.ab$loglik)
  
  # Allow a bit more tolerance in the coefficients
  expect_equal(fit.logbin$coefficients, fit.logbin.em$coefficients,
               tolerance = 1e-4)
  expect_equal(fit.logbin$coefficients, fit.logbin.glm$coefficients,
               tolerance = 1e-4)
  expect_equal(fit.logbin$coefficients, fit.logbin.glm2$coefficients,
               tolerance = 1e-4)
  expect_equal(fit.logbin$coefficients, fit.logbin.ab$coefficients,
               tolerance = 1e-4)
  
})

test_that("monotonicity constraints work", {
  
  require(glm2)
  data(heart)
  
  # Reverse the data so that constraints will be active
  heart$AgeGroup <- 4 - heart$AgeGroup
  heart$Severity <- 4 - heart$Severity
  
  # Should get a warning about the boundary
  expect_warning(
    fit.mono <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup) +
                         factor(Severity), data = heart, mono = 1:2,
                       start = c(-1, rep(c(0.1, 0.2), 2)))
  )
  
  # Estimate on the boundary
  expect_true(fit.mono$boundary)
  # Coefficients should be increasing
  expect_true(all(c(diff(c(0, fit.mono$coefficients[2:3])),
                    diff(c(0, fit.mono$coefficients[4:5]))) >= 0))
  
  # Also try with method = "em"
  expect_warning(
    fit.mono.em <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup) +
                            factor(Severity), data = heart, mono = 1:2,
                          method = "em", start = c(-1, rep(c(0.1, 0.2), 2)))
  )
  
  # Estimate on the boundary
  expect_true(fit.mono.em$boundary)
  # Coefficients should be increasing
  expect_true(all(c(diff(c(0, fit.mono.em$coefficients[2:3])),
                    diff(c(0, fit.mono.em$coefficients[4:5]))) >= 0))
  
  # Other methods do not support monotonicity
  expect_error(
    fit.mono.glm <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup) +
                             factor(Severity), data = heart, mono = 1:2,
                           method = "glm", start = c(-1, rep(c(0.1, 0.2), 2)))
  )
  expect_error(
    fit.mono.glm2 <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup) +
                              factor(Severity), data = heart, mono = 1:2,
                            method = "glm2", start = c(-1, rep(c(0.1, 0.2), 2)))
  )
  expect_error(
    fit.mono.ab <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup) +
                            factor(Severity), data = heart, mono = 1:2,
                          method = "ab", start = c(-1, rep(c(0.1, 0.2), 2)))
  )
  
  
})

test_that("offsets work", {
  
  require(glm2)
  data(heart)
  
  # Simple model so that we can compare with glm/glm2...
  # ...which don't work with the more full model
  
  fit.logbin <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup), 
                       data = heart, start = rep(-0.5, 3))
  
  # Put the same offset on all observations
  same_offset <- -0.3
  os <- rep(same_offset, nrow(heart))
  
  fit.logbin.os <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup),
                          data = heart, start = rep(-0.5, 3), offset = os)
  
  # Intercept should differ by the offset
  expect_equal(fit.logbin$coefficients[1],
               fit.logbin.os$coefficients[1] + same_offset)
  # Other coefficients should be the same
  expect_equal(fit.logbin$coefficients[-1],
               fit.logbin.os$coefficients[-1])
  
  # Don't allow positive offsets (for now?)
  expect_error(
    fit.logbin.osp <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup),
                             data = heart, start = rep(-0.5, 3), offset = -os)
  )
  
})