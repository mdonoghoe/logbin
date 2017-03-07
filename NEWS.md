# logbin 2.0.2

* Changed `logbin.ab` so that `outer.iterations` can be set by `control$maxit` in the call to `logbin`, and `maxit` (inner iterations) is set by `control.method$maxit`.
* `logbin` now returns `contrasts`, `qr`, `R` and `effects` components, so associated `glm` S3 methods (e.g. `influence`, `plot`) will work.
* `weights` component returned is the IRLS working weights at the last iteration, for consistency with `glm` (even though these are not actually used except with `method = "glm"`).
* `plot.logbin.smooth` with `type = "diagnostics"` will plot the model diagnostics usually produced by `plot.lm`.
* `plot.logbin.smooth` will allow extra parameters (e.g. `col`, `ylim`) to be passed.

# logbin 2.0.1

* (Had to increment version number to satisfy CRAN check after correcting an error)

# logbin 2.0

* `method` option added: fit using glm, glm2, adaptive barrier, CEM or single EM algorithm with overparameterised model
* `accelerate` option added: use `turboEM` methods to speed up convergence
* (Minor change) Corrected the definition of `boundary` 

# logbin 1.2

* Minor documentation changes (e.g. updated references)

# logbin 1.1

* Very minor documentation changes

# logbin 1.0

* Added `na.action` argument to `logbin()` and `smooth.logbin()`