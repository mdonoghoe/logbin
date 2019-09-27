logbin.cem <- function(mt, mf, Y, offset, mono, start, control, accelerate,
                       control.method, warn) {
  control2 <- control
  control2$trace <- (control$trace > 1)
  
  reparam <- logbin.reparam(mt, mf, "cem", mono)
  if (!is.null(start)) {
    start.expand <- logbin.expand(start, reparam, "cem")
  }
  
  best.model <- NULL
  best.loglik <- -Inf
  best.param <- NULL
  allconv <- TRUE
  totaliter <- 0
  
  if (control$coeftrace) Amat <- np.coefhist <- coefhist <- list()
    
  if(length(reparam$Vmat) == 0) {
    if(control$trace > 0) cat("logbin parameterisation 1/1\n")
    X <- model.matrix(mt, mf)
    best.model <- nplbin(Y, X, offset, start, control2,
                         accelerate, control.accelerate = list(control.method))
    best.loglik <- best.model$loglik
    best.param <- 0
    allconv <- best.model$converged
    totaliter <- best.model$iter
    if (control$coeftrace) np.coefhist[[1]] <- best.model$coefhist
    if (control$trace > 0 & control$trace <= 1)
      cat("Deviance =", best.model$deviance, "Iterations -", best.model$iter, "\n")
  } else {
    n.refvecs <- do.call("c", reparam$nref)
    n.param <- prod(n.refvecs)
    if (n.param > .Machine$integer.max)
      stop(paste0("The number of parameter subspaces is larger than R can handle ",
                  "(see .Machine$integer.max). Make some terms monotonic or use ",
                  "method = \"em\""))
    if (!is.null(start)) startrefs <- do.call("c", start.expand$vstar.id)
    else startrefs <- rep(1L, length(n.refvecs))
    designs.it <- do.call("iproduct", 
                          mapply(subsporder, npar = n.refvecs, 
                                 start = startrefs, SIMPLIFY = FALSE))
    paramcount <- 0L
    while (TRUE) {
      param <- try(iterators::nextElem(designs.it), silent = TRUE)
      if (inherits(param, "try-error") && param == "Error : StopIteration\n") break
      paramcount <- paramcount + 1L
      if (control$trace > 0) cat("logbin parameterisation ", paramcount, "/",
                                 n.param, "\n", sep = "")
      des <- logbin.design2(mt, mf, "cem", reparam, unlist(param))
      X <- des$X.reparam
      thismodel <- nplbin(Y, X, offset, if(!is.null(start) && paramcount == 1L) start.expand$coefs.exp else NULL,
                          control2, accelerate, control.accelerate = list(control.method))
      if (!thismodel$converged) allconv <- FALSE
      if (control$coeftrace) {
        np.coefhist[[paramcount]] <- thismodel$coefhist
        Amat[[paramcount]] <- des$A
      }
      if (control$trace > 0 & control$trace <= 1)
        cat("Deviance =", thismodel$deviance, "Iterations -", thismodel$iter, "\n")
      totaliter <- totaliter + thismodel$iter
      if(thismodel$loglik > best.loglik) {
        best.model <- thismodel
        best.loglik <- thismodel$loglik
        best.param <- unlist(param)
        if(thismodel$converged & !thismodel$boundary) break
      }
    }
    
    #n.refvecs <- do.call("c", reparam$nref)
    #design.all <- expand.grid(lapply(n.refvecs, seq_len))
    #nparam <- nrow(design.all)
    
    #for(param in seq_len(nparam)) {
    #   if(control$trace > 0) cat("logbin parameterisation ",param,"/",nparam,"\n",sep="")
    #   des <- logbin.design2(mt, mf, "cem", reparam, design.all[param,])
    #   X <- des$X.reparam
    #   thismodel <- nplbin(Y, X, offset, NULL,
    #                       control2, accelerate, control.accelerate = list(control.method))
    #   if (!thismodel$converged) allconv <- FALSE
    #   if (control$coeftrace) np.coefhist[[param]] <- thismodel$coefhist
    #   if (control$trace > 0 & control$trace <= 1)
    #     cat("Deviance =", thismodel$deviance, "Iterations -", thismodel$iter, "\n")
    #   totaliter <- totaliter + thismodel$iter
    #   if(thismodel$loglik > best.loglik) {
    #     best.model <- thismodel
    #     best.loglik <- thismodel$loglik
    #     best.param <- param
    #     if(thismodel$converged & !thismodel$boundary) break
    #   }
    # }
  }
    
  if (length(reparam$Vmat) == 0) {
    np.coefs <- coefs <- coefs.boundary <- best.model$coefficients
    nn.design <- design <- X
    if (control$coeftrace) coefhist[[1]] <- np.coefhist
  } else {
    np.coefs <- best.model$coefficients
    best.design <- logbin.design2(mt, mf, "cem", reparam, best.param)
    nn.design <- best.design$X.reparam
    coefs <- as.vector(logbin.reduce(np.coefs, best.design$A))
    names(coefs) <- gsub("`", "", colnames(best.design$X.orig))
    design <- best.design$X.orig
    coefs.boundary <- np.coefs[logbin.expand(coefs, reparam, "cem")$which.boundary]
    if (control$coeftrace) {
      coefhist <- mapply(coefhist.reduce, np.coefhist = np.coefhist, Amat = Amat, 
                         MoreArgs = list(cnames = names(coefs)),
                         SIMPLIFY = FALSE)
    }
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
              fitted.values = best.model$fitted.values,
              linear.predictors = best.model$linear.predictors, deviance = best.model$deviance,
              loglik = best.model$loglik, aic = best.model$aic, aic.c = best.model$aic.c,
              null.deviance = best.model$null.deviance, iter = c(totaliter, best.model$iter),
              prior.weights = best.model$prior.weights,
              df.residual = best.model$df.residual, df.null = best.model$df.null,
              y = best.model$y, x = design, converged = best.model$converged,
              boundary = boundary, np.coefficients = np.coefs,
              nn.x = nn.design)
  if (control$coeftrace) fit$coefhist <- coefhist
  fit
}

# This is a function that determines the ordering of the parameter subspaces
subsporder <- function(npar, start) {
  ((seq_len(npar) + start - 2L) %% npar) + 1L
}