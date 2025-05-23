# logbin 2.0.6.9000 (under development)

...

# logbin 2.0.6

* Replace `vcov` calculation in `summary.logbin` to replicate the version in `summary.glm`.
* Implement tests using `testthat`.
* Clarify the copyright status of functions based on base R.
* (Minor) Removed unused `contr.isototonic.rev` function. Also no longer exporting `interpret.logbin.smooth`.
* (Minor) Suppress startup messages when `turboEM::turboem` is called for the first time.
* (Minor) Update `logbin-package` Rd file to remove CRAN `NOTE`.

# logbin 2.0.5

* Fixed reparameterisation for factor variables in `method = "cem"` (@jakobschoepe, #3)
* Implemented a more efficient approach to parameter expansion (using the transformation matrix).
* When using the parameter expansion approach, check convergence of the reduced parameter vector. As an indirect result, the `conv.test` function is no longer exported.
* Fixed error in `logbin` if offset provided
* Implement `coeftrace` via `turboEM::turboem` argument `keep.paramval`, to avoid taking a different path when it is turned on (#5)

# logbin 2.0.4

* Added citation to _Journal of Statistical Software_ article.

# logbin 2.0.3

* Added `coeftrace` argument to `logbin.control`, such that coefficient history can be obtained for `method = "em"` and `method = "cem"`.

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
* (Minor) Corrected the definition of `boundary` 

# logbin 1.2

* Minor documentation changes (e.g. updated references)

# logbin 1.1

* Very minor documentation changes

# logbin 1.0

* Added `na.action` argument to `logbin()` and `smooth.logbin()`