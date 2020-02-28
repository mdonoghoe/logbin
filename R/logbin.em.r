logbin.em <- function(mt, mf, Y, offset, mono, start, control, accelerate = c("em","squarem","pem","qn"),
                      control.method, warn)
{
  accelerate = match.arg(accelerate)
  control2 <- control
  control2$trace <- (control$trace > 1)
  
  reparam <- logbin.reparam(mt, mf, "em", mono)
  if (!is.null(start)) {
    start.expand <- logbin.expand(start, reparam, "em")
  }
  
  if(control$trace > 0) cat("logbin parameterisation 1/1\n")
  if (length(reparam$Vmat) == 0) {
    X <- model.matrix(mt, mf)
    Amat <- diag(ncol(X))
  }
  else {
    des <- logbin.design2(mt, mf, "em", reparam)
    X <- des$X.reparam
    Amat <- des$A
  }
  
  thismodel <- nplbin(Y, X, offset, if (!is.null(start)) start.expand$coefs.exp else NULL, 
                      Amat = Amat, control = control2, accelerate = accelerate, 
                      control.accelerate = list(control.method))
  
  if(control$trace > 0 & control$trace <= 1)
    cat("Deviance =", thismodel$deviance, "Iterations -", thismodel$iter, "\n")
  
  if (length(reparam$Vmat) == 0) {
    np.coefs <- coefs <- coefs.boundary <- thismodel$coefficients
    nn.design <- design <- X
    if (control$coeftrace) coefhist <- thismodel$coefhist
  } else {
    np.coefs <- thismodel$coefficients
    nn.design <- X
    #reparam <- logbin.reparameterise(np.coefs, mt, mf, "em", allref$allref, allref$monotonic, design.all[1,])
    coefs <- as.vector(logbin.reduce(np.coefs, des$A))
    names(coefs) <- gsub("`", "", colnames(des$X.orig))
    design <- des$X.orig
    coefs.boundary <- np.coefs[logbin.expand(coefs, reparam, "em")$which.boundary]
    if (control$coeftrace) {
      coefhist <- coefhist.reduce(thismodel$coefhist, des$A, names(coefs))
    }
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
              fitted.values = thismodel$fitted.values,
              linear.predictors = thismodel$linear.predictors, deviance = thismodel$deviance,
              loglik = thismodel$loglik, aic = thismodel$aic - 2*vardiff, aic.c = aic.c,
              null.deviance = thismodel$null.deviance, iter = thismodel$iter,
              prior.weights = thismodel$prior.weights,
              df.residual = thismodel$df.residual + vardiff, df.null = thismodel$df.null,
              y = thismodel$y, x = design, converged = thismodel$converged,
              boundary = boundary, np.coefficients = np.coefs,
              nn.x = nn.design)
  if (control$coeftrace) fit$coefhist <- coefhist
  fit
}