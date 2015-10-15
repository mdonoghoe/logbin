utils::globalVariables("tol")

nplbin <- function(y, x, offset, start, control = logbin.control(), 
                   accelerate = c("em","squarem","pem","qn"), control.accelerate = list(list())) {
  control <- do.call("logbin.control", control)
  accelerate <- match.arg(accelerate)
  
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y))
    rownames(y)
  else names(y)
  
  if (any(x < 0)) stop("x must be non-negative")
  if (any(apply(x, 2, function(col) all(col==0)))) stop("x contains column with all 0")
  
  nvars <- ncol(x)
  nobs <- NROW(y)
  
  n <- weights <- rep(1, nobs)
  if (is.null(offset)) offset <- rep.int(0, nobs)
  if (any(offset > 0))
    stop("offset must be non-positive")
  
  fam <- binomial(link = log)
  eval(fam$initialize)
  
  mu.eta <- fam$mu.eta
  linkinv <- fam$linkinv
  dev.resids <- fam$dev.resids
  aic <- fam$aic
    
  y1 <- round(n*y)
  y2 <- round(n*(1-y))
  
  x.mins <- apply(x, 2, function(t) min(t[t > 0]))
  x.scale <- 1/x.mins
  
  converged <- FALSE
    
  coefold <- if (!is.null(start)) {
    if (length(start) != nvars)
      stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                    nvars, paste(deparse(xnames), collapse = ", ")), domain = NA)
    else if (any(start >= -control$bound.tol))
      stop("'start' is on our outside the boundary of the parameter space (consider 'bound.tol')", domain = NA)
    else start
  } else {
    simple <- log(mean(y)) / colMeans(x) - 2 * control$bound.tol
    logy <- log(y)
    logy[is.infinite(logy)] <- min(logy[!is.infinite(logy)]) - 2 * control$bound.tol
    trymat <- tryCatch(as.vector(solve(t(x) %*% x) %*% t(x) %*% (logy)) + 2*control$bound.tol,
                       error = function(e) NULL)
    if (is.null(trymat)) simple
    else if (any(trymat >= -control$bound.tol)) simple
    else trymat
  }
    
  fixptfn <- function(p, y1, y2, n, x, x.s, o, nobs, nvars, fam, bound.tol) {
    x1 <- sweep(x, 2, FUN = "*", x.s)
    p.scale <- p / x.s
    eta <- drop(x1 %*% p.scale) + o
    estep <- y1 + y2 * ((matrix(fam$linkinv(p.scale), nobs, nvars, byrow = TRUE) - fam$linkinv(eta))/(1 - fam$linkinv(eta)))
    pnew.scale <- log(colSums(estep * x1) / colSums(n * x1))
    pnew <- pnew.scale * x.s
    pnew[pnew >= 0] <- -bound.tol / 2
    return(pnew)
  }
  
  objfn <- function(p, y1, y2, n, x, x.s, o, nobs, nvars, fam, bound.tol) {
    eta <- drop(x %*% p) + o
    mu <- n * fam$linkinv(eta)
    negll <- -sum(dbinom(y1, size = n, prob = mu / n, log = TRUE))
    return(negll)
  }
  
  validparams <- function(p) return(all(p <= 0))
  
  conv.user <- function(old,new) return(conv.test(old, new, tol))
  
  res <- turboEM::turboem(par = coefold, fixptfn = fixptfn, objfn = objfn, method = accelerate,
                          pconstr = validparams, y1 = y1, y2 = y2, n = n, x = x, x.s = x.scale,
                          o = offset, nobs = nobs, nvars = nvars, fam = fam, bound.tol = control$bound.tol,
                          control.run = list(convtype = "parameter", tol = control$epsilon,
                                             stoptype = "maxiter", maxiter = control$maxit,
                                             convfn.user = conv.user, trace = control$trace),
                          control.method = control.accelerate)
  coefnew <- res$pars[1,]
  names(coefnew) <- xnames
    
  eta <- drop(x %*% coefnew) + offset
  mu <- n * linkinv(eta)
  residuals <- (y - (mu / n)) / mu.eta(eta)
    
  names(y) <- names(mu) <- names(eta) <- names(residuals) <- ynames
    
  dev.new <- sum(dev.resids(y, mu / n, n))
  aic.model <- aic(y, n, mu / n, weights, dev.new) + 2 * nvars
  aic.c.model <- aic.model + 2 * nvars * (nvars + 1) / (nobs - nvars - 1)
    
  wtdmu <- sum(n * y) / sum(n)
  nulldev <- sum(dev.resids(y, wtdmu, n))
  nulldf <- nobs - 1
  resdf <- nobs - nvars
    
  boundary <- any(coefnew > -control$bound.tol)
    
  list(coefficients = coefnew, residuals = residuals, fitted.values = mu / n, rank = nvars,
       family = fam, linear.predictors = eta, deviance = dev.new, aic = aic.model, 
       aic.c = aic.c.model, null.deviance = nulldev, iter = res$itr[1], weights = weights, 
       prior.weights = n, df.residual = resdf, df.null = nulldf, y = y, 
       converged = res$convergence[1], boundary = boundary, loglik = -res$value.objfn[1], 
       nn.design = x)   
}