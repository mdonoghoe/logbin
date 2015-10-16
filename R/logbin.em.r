logbin.em <- function(mt, mf, Y, offset, mono, start, control, accelerate = c("em","squarem","pem","qn"),
                      control.method, warn)
{
  accelerate = match.arg(accelerate)
  control2 <- control
  control2$trace <- (control$trace > 1)

  allref <- logbin.allref(mt, mf, "em", mono, start)
  design.numref <- sapply(allref$allref, length)
  
  if (length(allref$allref) == 0)
    X <- model.matrix(allref$terms, allref$data)
  else {
    design.all <- expand.grid(lapply(design.numref, seq_len))
    X <- logbin.design(allref$terms, allref$data, "em", allref$allref, allref$monotonic, design.all[1,])
  }
  
  thismodel <- nplbin(Y, X, offset, allref$start.new, control2, accelerate, control.accelerate = list(control.method))
  
  if (length(allref$allref) == 0) {
    np.coefs <- coefs <- coefs.boundary <- thismodel$coefficients
    nn.design <- design <- X
  } else {
    np.coefs <- thismodel$coefficients
    nn.design <- X
    reparam <- logbin.reparameterise(np.coefs, mt, mf, "em", allref$allref, allref$monotonic, design.all[1,])
    coefs <- reparam$coefs
    design <- reparam$design
    coefs.boundary <- reparam$coefs.boundary
  }
  
  nvars <- length(coefs)
  vardiff <- length(np.coefs) - nvars
  aic.c <- thismodel$aic - 2 * vardiff + 2 * nvars * (nvars + 1) / (NROW(Y) - nvars - 1)
  
  boundary <- any(coefs.boundary > -control$bound.tol)
  
  if (warn) {
    if (!thismodel$converged) {
      if (identical(accelerate, "em"))
        warning(gettextf("nplbin: algorithm did not converge within %d iterations -- increase 'maxit'.",
                         control$maxit), call. = FALSE)
      else
        warning(gettextf("nplbin(%s): algorithm did not converge within %d iterations -- increase 'maxit' or try with 'accelerate = \"em\"'.",
                         accelerate, control$maxit), call. = FALSE)
    }
    if(boundary) {
      if(coefs.boundary[1] > -control$bound.tol)
          warning("nplbin: fitted probabilities numerically 1 occurred", call. = FALSE)
      else
          warning("nplbin: MLE on boundary of constrained parameter space", call. = FALSE)
    }
  }
  
  fit <- list(coefficients = coefs, residuals = thismodel$residuals,
              fitted.values = thismodel$fitted.values, rank = nvars,
              linear.predictors = thismodel$linear.predictors, deviance = thismodel$deviance,
              loglik = thismodel$loglik, aic = thismodel$aic - 2*vardiff, aic.c = aic.c,
              null.deviance = thismodel$null.deviance, iter = thismodel$iter,
              prior.weights = thismodel$prior.weights, weights = rep(1, NROW(Y)),
              df.residual = thismodel$df.residual + vardiff, df.null = thismodel$df.null,
              y = thismodel$y, x = design, converged = thismodel$converged,
              boundary = boundary, np.coefficients = np.coefs,
              nn.x = nn.design)
  fit
}