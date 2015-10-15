logbin.ab <- function(mt, mf, Y, offset, mono, start, control, control.method, warn) {
  x <- as.matrix(model.matrix(mt, mf))
  xnames <- dimnames(x)
  ynames <- if (is.matrix(y))
    rownames(y)
  else names(y)
  
  nvars <- ncol(X)
  nobs <- NROW(y)
  
  n <- rep(1, nobs)
  if (is.null(offset)) offset <- rep.int(0, nobs)
  
  fam <- binomial(link = log)
  eval(fam$initialize)
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
  
  negll <- function(theta, y, n, x, o) {
    p.fit <- exp(drop(x %*% theta) + o)
    -sum(dbinom(y, size = n, prob = p.fit, log = TRUE))
  }
  
  gradll <- function(theta, y, n, x, o) {
    # Gradient of **negative** log-likelihood
  }
  
  fit.ab <- constrOptim(theta.start, f = negll, y = y1, n = n, x = x, o = offset)
}