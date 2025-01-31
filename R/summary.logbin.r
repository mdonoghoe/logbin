#  Based on File src/library/stats/R/glm.R
#  Part of the R package, https://www.R-project.org
#
#  Modified by Mark W. Donoghoe
#     31/01/2025 - to work with logbin objects
#
#  Copyright (C) 1995-2020 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

summary.logbin <- function(object, correlation = FALSE, ...) {
  df.r <- object$df.residual
  p <- object$rank
  coef.p <- object$coefficients
  
  dispersion <- 1
  aliased <- is.na(coef(object))
  
  if(!object$boundary) {
    x <- object$x
    s <- object$prior.weights
    y <- s * object$y
    mu <- object$fitted.values
    if (p > 0) {
      p1 <- 1L:p
      Qr <- object$qr
      coef.p <- object$coefficients[Qr$pivot[p1]]
      covmat.unscaled <- chol2inv(Qr$qr[p1,p1,drop=FALSE])
      dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
      covmat <- dispersion*covmat.unscaled
      var.cf <- diag(covmat)
      s.err <- sqrt(var.cf)
      tvalue <- coef.p/s.err
      pvalue <- 2 * pnorm(-abs(tvalue))
      coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
      if (!is.null(object$call$mono))
        warning("model contains monotonicity constraints, asymptotic covariance matrix may not be valid", call. = FALSE)
    } else {
      coef.table <- matrix(, 0L, 4L)
      covmat.unscaled <- covmat <- matrix(, 0L, 0L)
      df.f <- length(aliased)
    }
  }
  else {
    warning("MLE on boundary of parameter space, cannot use asymptotic covariance matrix", call. = FALSE)
    covmat.unscaled <- matrix(NaN, p, p)
    covmat <- matrix(NaN, p, p)
    dimnames(covmat.unscaled) <- dimnames(covmat) <- list(names(coef.p), names(coef.p))
    coef.table <- cbind(coef.p, NaN, NaN, NaN)
  }
  
  dimnames(coef.table) <- list(names(coef.p), c("Estimate","Std. Error","z value","Pr(>|z|)"))
  
  keep <- match(c("call", "family", "deviance", "aic", "aic.c", "df.residual",
                  "null.deviance", "df.null", "iter", "na.action", "method"), names(object), 0L)
  ans <- c(object[keep], list(deviance.resid = residuals(object,type="deviance"),
                              coefficients = coef.table, aliased = aliased,
                              dispersion = dispersion, df = c(p, df.r, p),
                              cov.unscaled = covmat.unscaled, cov.scaled = covmat))
  if(correlation && !any(is.nan(covmat.unscaled))) {
    dd <- sqrt(diag(covmat.unscaled))
    ans$correlation <- covmat.unscaled/outer(dd, dd)
  }
  if(inherits(object,"logbin.smooth")) ans$knots <- object$knots
  class(ans) <- c("summary.logbin", "summary.glm")
  ans
}