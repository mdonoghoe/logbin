logbin.cem <- function(mt, mf, Y, offset, mono, start, control, accelerate,
                       control.method, warn) {
  control2 <- control
  control2$trace <- (control$trace > 1)

  allref <- logbin.allref(mt, mf, "cem", mono, start)
  design.numref <- sapply(allref$allref, length)
  
  best.model <- NULL
  best.loglik <- -Inf
  best.param <- NULL
  allconv <- TRUE
  totaliter <- 0
    
  if(length(allref$allref) == 0) {
    if(control$trace > 0) cat("logbin parameterisation 1/1\n")
    X <- model.matrix(allref$terms, allref$data)
    best.model <- nplbin(Y, X, offset, allref$start.new, control2,
                         accelerate, control.accelerate = list(control.method))
    best.loglik <- best.model$loglik
    best.param <- 0
    allconv <- best.model$converged
    totaliter <- best.model$iter
    if(control$trace > 0 & control$trace <= 1)
      cat("Deviance =", best.model$deviance, "Iterations -", best.model$iter, "\n")
  } else {
    design.all <- expand.grid(lapply(design.numref, seq_len))
    nparam <- nrow(design.all)
    
    for(param in seq_len(nparam)) {
      if(control$trace > 0) cat("logbin parameterisation ",param,"/",nparam,"\n",sep="")
      X <- logbin.design(allref$terms, allref$data, "cem", allref$allref, allref$monotonic, design.all[param,])
      thismodel <- nplbin(Y, X, offset, if (param == 1) allref$start.new else NULL,
                          control2, accelerate, control.accelerate = list(control.method))
      if(!thismodel$converged) allconv <- FALSE
      if(control$trace > 0 & control$trace <= 1)
        cat("Deviance =", thismodel$deviance, "Iterations -", thismodel$iter, "\n")
      totaliter <- totaliter + thismodel$iter
      if(thismodel$loglik > best.loglik) {
        best.model <- thismodel
        best.loglik <- thismodel$loglik
        best.param <- param
        if(thismodel$converged & !thismodel$boundary) break
      }
    }
  }
    
  if (length(allref$allref) == 0) {
    np.coefs <- coefs <- coefs.boundary <- best.model$coefficients
    nn.design <- design <- model.matrix(allref$terms, allref$data)
  } else {
    np.coefs <- best.model$coefficients
    nn.design <- logbin.design(allref$terms, allref$data, "cem", allref$allref, allref$monotonic, design.all[best.param,])
    reparam <- logbin.reparameterise(np.coefs, mt, mf, "cem", allref$allref, allref$monotonic, design.all[best.param,])
    coefs <- reparam$coefs
    design <- reparam$design
    coefs.boundary <- reparam$coefs.boundary
  }
  
  boundary <- any(coefs.boundary > -control$bound.tol)
  
  if (warn) {
    if (!best.model$converged | (!allconv & best.model$boundary)) {
      if (identical(accelerate, "em"))
        warning(gettextf("nplbin: algorithm did not converge within %d iterations -- increase 'maxit'.", control$maxit), 
                call. = FALSE)
      else
        warning(gettextf("nplbin(%s): algorithm did not converge within %d iterations -- increase 'maxit' or try with 'accelerate = \"em\"'.",
                         accelerate, control$maxit),
                call. = FALSE)
    }
    if (boundary) {
      if (coefs.boundary[1] > -control$bound.tol)
          warning("nplbin: fitted probabilities numerically 1 occurred", call. = FALSE)
      else
          warning("nplbin: MLE on boundary of constrained parameter space", call. = FALSE)
    }
  }
    
  fit <- list(coefficients = coefs, residuals = best.model$residuals, 
              fitted.values = best.model$fitted.values, rank = best.model$rank,
              linear.predictors = best.model$linear.predictors, deviance = best.model$deviance,
              loglik = best.model$loglik, aic = best.model$aic, aic.c = best.model$aic.c,
              null.deviance = best.model$null.deviance, iter = c(totaliter, best.model$iter),
              prior.weights = best.model$prior.weights, weights = rep(1, NROW(Y)),
              df.residual = best.model$df.residual, df.null = best.model$df.null,
              y = best.model$y, x = design, converged = best.model$converged,
              boundary = boundary, np.coefficients = np.coefs,
              nn.x = nn.design)
  fit
}