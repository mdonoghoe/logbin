
<!-- README.md is generated from README.Rmd. Please edit that file -->
logbin
======

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/logbin)](http://cran.r-project.org/package=logbin) [![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/logbin)](http://cran.rstudio.com/web/packages/logbin/index.html)

`logbin` provides methods for performing relative risk regression by fitting log-link GLMs and GAMs to binomial data. As well as providing a consistent interface to use the usual Fisher scoring algorithm (via `glm` or `glm2`) and an adaptive barrier approach (via `constrOptim`), it implements EM-type algorithms that have more stable convergence properties than other methods.

An example of periodic non-convergence using `glm` (run with `trace = TRUE` to see deviance at each iteration):

``` r
require(glm2, quietly = TRUE)
data(heart)

start.p <- sum(heart$Deaths) / sum(heart$Patients)
t.glm <- system.time(
  fit.glm <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup) + factor(Severity) + 
                      factor(Delay) + factor(Region), data = heart,
                    start = c(log(start.p), -rep(1e-4, 8)), method = "glm", maxit = 10000)
)
```

The combinatorial EM method (Marschner and Gillett, 2012) provides stable convergence:

``` r
t.cem <- system.time(fit.cem <- update(fit.glm, method = "cem"))
```

...but it can take a while. Using an overparameterised EM approach removes the need to run \(3^4 = 81\) separate EM algorithms:

``` r
t.em <- system.time(fit.em <- update(fit.glm, method = "em"))
```

...while generic EM acceleration algorithms (from the `turboEM` package) can speed this up further still:

``` r
t.cem.acc <- system.time(fit.cem.acc <- update(fit.cem, accelerate = "squarem"))
t.em.acc <- system.time(fit.em.acc <- update(fit.em, accelerate = "squarem"))
```

Comparison of results:

    #>         converged    logLik iterations   time
    #> glm         FALSE -186.7366      10000   2.22
    #> cem          TRUE -179.9016     445161 109.29
    #> em           TRUE -179.9016       7403   1.91
    #> cem.acc      TRUE -179.9016       7974   9.21
    #> em.acc       TRUE -179.9016         92   0.13

<!-- Still to add: adaptive barrier, semi-parametric examples -->
References
----------

-   Marschner, I. C. and A. C. Gillett (2012). Relative risk regression: reliable and flexible methods for log-binomial models. *Biostatistics* **13**(1): 179–192.
