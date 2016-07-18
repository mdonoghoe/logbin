
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
    #> glm         FALSE -186.7366      10000   2.21
    #> cem          TRUE -179.9016     445161 103.13
    #> em           TRUE -179.9016       7403   1.84
    #> cem.acc      TRUE -179.9016       7974   8.72
    #> em.acc       TRUE -179.9016         92   0.10

An adaptive barrier algorithm can also be applied using `method = "ab"`, with user-specified options via `control.method`: see `help(logbin)` for more details.

Semi-parametric regression using B-splines (Donoghoe and Marschner, 2015) can be incorporated by using the `logbin.smooth` function. See `example(logbin.smooth)` for a simple example.

Installation
------------

Get the released version from CRAN:

``` r
install.packages("logbin")
```

Or the development version from github:

``` r
# install.packages("devtools")
devtools::install_github("mdonoghoe/logbin")
```

References
----------

-   Donoghoe, M. W. and I. C. Marschner (2015). Flexible regression models for rate differences, risk differences and relative risks. *International Journal of Biostatistics* **11**(1): 91-108.
-   Marschner, I. C. and A. C. Gillett (2012). Relative risk regression: reliable and flexible methods for log-binomial models. *Biostatistics* **13**(1): 179-192.
