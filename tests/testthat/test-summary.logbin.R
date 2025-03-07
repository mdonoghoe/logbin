test_that("summary.logbin gives same results as summary.glm", {
  
  require(glm2)
  data(heart)
  fit.glm <- glm(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup),
                 family = binomial(log), data = heart,
                 start = rep(-0.5, 3))
  
  fit.logbin <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup), 
                       data = heart, start = rep(-0.5, 3), method = "glm")
  
  s.glm <- summary(fit.glm)
  s.logbin <- summary(fit.logbin)
  
  common_parts <- c("family", "deviance", "aic", "df.residual", "null.deviance",
                    "df.null", "iter", "deviance.resid", "coefficients",
                    "aliased", "dispersion", "df", "cov.unscaled", "cov.scaled")
  
  expect_equal(s.glm[common_parts], s.logbin[common_parts])
  
  expect_equal(vcov(fit.glm), vcov(fit.logbin))
  expect_equal(confint.default(fit.glm), confint(fit.logbin))
  
})

test_that("NaN covariance matrix for boundary estimate", {
  
  dat.boundary <- data.frame(y = 10, n = 10)
  # Converges to boundary estimate
  expect_warning(mdl <- logbin(cbind(y, n-y) ~ 1, data = dat.boundary))
  # Gives NaN covariance matrix & confidence intervals
  nan_vcov <- matrix(NaN, 1, 1, dimnames = list("(Intercept)", "(Intercept)"))
  nan_confint <- matrix(NaN, 1, 2, dimnames = list("(Intercept)", c("2.5 %", "97.5 %")))
  expect_warning(v <- vcov(mdl))
  expect_identical(v, nan_vcov)
  expect_warning(ci <- confint(mdl))
  expect_identical(ci, nan_confint)
  
})