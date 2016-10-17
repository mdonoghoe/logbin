logbin.ab <- function(mt, mf, Y, offset, mono, start, control, control.method, warn) {
  x <- as.matrix(model.matrix(mt, mf))
  xnames <- dimnames(x)[[2L]]
  y <- Y
  ynames <- if (is.matrix(y))
    rownames(y)
  else names(y)
  
  nvars <- ncol(x)
  nobs <- NROW(y)
  
  n <- weights <- rep(1, nobs)
  if (is.null(offset)) offset <- rep.int(0, nobs)
  
  fam <- binomial(link = log)
  eval(fam$initialize)
  
  mu.eta <- fam$mu.eta
  linkinv <- fam$linkinv
  dev.resids <- fam$dev.resids
  aic <- fam$aic
  
  y1 <- round(n * y)
  
  if (!is.null(start)) {
    if (length(start) != nvars)
      stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                    nvars, paste(deparse(xnames), collapse = ", ")),
           domain = NA)
    else if (any(x %*% start + offset >= 0))
      stop("'start' is on or outside the boundary of the parameter space", domain = NA)
    else
      theta.start <- start
  } else {
    allref <- logbin.allref(mt, mf, "cem", mono, NULL)
    if (length(allref$allref) == 0)
      theta.start <- log(mean(y)) - 2 * control$bound.tol
    else {
      design.numref <- sapply(allref$allref, length)
      design.all <- expand.grid(lapply(design.numref, seq_len))
      start.np <- rep((log(mean(y)) - 2 * control$bound.tol) / nvars, nvars)
      reparam <- logbin.reparameterise(start.np, mt, mf, "cem", allref$allref, allref$monotonic, design.all[1, ])
      theta.start <- reparam$coefs
    }
  }
  
  negll <- function(theta, y, n, x, offset) {
    p.fit <- exp(drop(x %*% theta) + offset)
    -sum(dbinom(y, size = n, prob = p.fit, log = TRUE))
  }
  
  gradll <- function(theta, y, n, x, offset) {
    p.fit <- exp(drop(x %*% theta) + offset)
    ll.grad <- y * x - (n - y) * x * (p.fit / (1 - p.fit))
    ll.grad[p.fit %in% c(0, 1), ] <- 0
    -colSums(ll.grad)
  }
  
  if (!is.null(control.method$method)) {
    method <- control.method$method
    control.method$method <- NULL
  } else
    method <- "BFGS"
  
  control.method$trace <- pmax(control$trace - 1, 0)
  
  fit.ab <- constrOptim(theta.start, f = negll, grad = gradll, ui = -x, ci = 0, 
                        control = control.method, method = method, outer.iterations = control$maxit,
                        outer.eps = control$epsilon, y = y1, n = n, x = x, offset = offset,
                        hessian = FALSE)
  
  coefficients <- fit.ab$par
  names(coefficients) <- xnames
  eta <- drop(x %*% coefficients) + offset
  mu <- n * linkinv(eta)
  residuals <- (y - (mu / n)) / mu.eta(eta)
  deviance <- sum(dev.resids(y, mu / n, n))
  loglik <- -fit.ab$value
  aic.model <- aic(y, n, mu / n, weights, dev.new) + 2 * nvars
  aic.c <- aic.model + 2 * nvars * (nvars + 1) / (nobs - nvars - 1)
  wtdmu <- sum(n * y) / sum(n)
  null.deviance <- sum(dev.resids(y, wtdmu, n))
  iter <- fit.ab$counts["function"]
  names(iter) <- NULL
  
  converged <- as.logical(fit.ab$convergence == 0)
  boundary <- any(eta >= -control$bound.tol)
  
  if (warn) {
    if (!converged) {
      if (fit.ab$convergence == 10)
        warning("constrOptim: algorithm did not converge -- Nelder-Mead simplex degenerate", call. = FALSE)
      else
        warning(paste0("constrOptim: algorithm did not converge (message: ",
                       fit.ab$message, ")"), call. = FALSE)
    }
    if (boundary)
      warning("constrOptim: fitted probabilities numerically 1 occurred", call. = FALSE)
  }
  
  list(coefficients = coefficients, residuals = residuals, fitted.values = mu / n,
       rank = nvars, linear.predictors = eta, deviance = deviance,
       loglik = -fit.ab$value, aic = aic.model, aic.c = aic.c, null.deviance = null.deviance,
       iter = iter, prior.weights = n,
       df.residual = nobs - nvars, df.null = nobs - 1, converged = converged,
       boundary = boundary)
}