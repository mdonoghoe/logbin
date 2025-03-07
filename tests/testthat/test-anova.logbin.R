test_that("anova.logbin doesn't work for single model", {
  
  require(glm2)
  data(heart)
  fit.logbin <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup),
                       data = heart,
                       start = rep(-0.5, 3))
  
  expect_error(anova(fit.logbin))
  
})

test_that("anova.logbin gives same result as anova.glm", {
  
  require(glm2)
  data(heart)
  fit.glm <- glm(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup),
                 family = binomial(log), data = heart,
                 start = rep(-0.5, 3))
  fit.glm2 <- glm(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup) + factor(Region),
                  family = binomial(log), data = heart,
                  start = rep(-0.5, 5))
  
  fit.logbin <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup), 
                       data = heart, start = rep(-0.5, 3), method = "glm")
  fit.logbin2 <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup) + factor(Region), 
                        data = heart, start = rep(-0.5, 5), method = "glm")
  
  expect_identical(anova(fit.glm, fit.glm2, test = "LRT"),
                   anova(fit.logbin, fit.logbin2, test = "LRT"))
  
})