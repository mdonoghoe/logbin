logbin <- function (formula, mono = NULL, data, subset, na.action, start = NULL, offset,
            control = list(...), model = TRUE, method = c("cem", "em", "glm", "glm2", "ab"),
            accelerate = c("em", "squarem", "pem", "qn"), control.method = list(), 
            warn = TRUE, ...) {
  call <- match.call()
  method <- match.arg(method)
  accelerate <- match.arg(accelerate)
  
  family <- binomial(link = log)
  if(missing(data)) data <- environment(formula)
  control <- do.call("logbin.control", control)
    
  outnames <- c("coefficients", "residuals", "fitted.values", "effects", "R", "rank", "qr", "family",
                "linear.predictors", "deviance", "loglik", "aic", "aic.c",
                "null.deviance", "iter", "prior.weights", "weights",
                "df.residual", "df.null", "y", "x")
  if (model) outnames <- c(outnames, "model")
  outnames <- c(outnames, "converged", "boundary", "na.action", "call",
                "formula", "terms", "data", "offset", "control", "method", "contrasts", "xlevels",
                "xminmax", "np.coefficients", "nn.x")
  if (control$coeftrace) outnames <- c(outnames, "coefhist")
  fit <- sapply(outnames, function(x) NULL)
  
  if (method %in% c("glm", "glm2")) {
    if (!is.null(mono)) stop(paste(method, "does not support monotonic constraints"))
    if (!identical(accelerate, "em")) warning(paste(method, "does not support acceleration"))
    res <- do.call(method, list(formula = call$formula, family = family, data = data,
                                start = call$start, offset = call$offset, 
                                subset = call$subset,
                                control = glm.control(epsilon = control$epsilon, 
                                                      maxit = control$maxit, 
                                                      trace = (control$trace > 1)), 
                                model = model, x = TRUE, y = TRUE, contrasts = NULL))
    ll.res <- logLik(res)
    mres <- match(outnames, names(res), 0L)
    fit[names(res)[mres]] <- res[mres]
    fit$loglik <- as.numeric(ll.res)
    fit$aic.c <- AIC(res, k = 2 * nobs(res) / (nobs(res) - attr(ll.res, "df") - 1))  
  } else {
    if (method == "ab") {
      if (!is.null(mono)) stop(paste(method, "does not currently support monotonic constraints"))
      if (!identical(accelerate, "em")) warning(paste(method, "does not support acceleration"))
    }
  
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    
    if (is.empty.model(mt)) stop("empty model")
    if (attr(mt, "intercept") != 1) stop("models without intercept are not supported by logbin")
    if (any(attr(mt, "order") > 1)) stop("models with interactions are not supported by logbin")
    if (attr(mt, "response") == 0) stop("missing response")
    
    offset <- as.vector(model.offset(mf))
    if(!is.null(offset)) {
      if(length(offset) != NROW(Y))
        stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                      length(offset), NROW(Y)), domain = NA)
    }
    
    Y <- model.response(mf, "numeric")
    if(length(dim(Y)) == 1L) {
      nm <- rownames(Y)
      dim(Y) <- NULL
      if(!is.null(nm)) names(Y) <- nm
    }
  
    logbin.method <- paste("logbin", method, sep = ".")
    logbin.args <- list(mt = mt, mf = mf, Y = Y, offset = offset, mono = mono,
                        start = start, control = control)
    if (method %in% c("cem", "em")) logbin.args$accelerate <- accelerate
    logbin.args <- c(logbin.args, list(control.method = control.method, warn = warn))
    res <- do.call(logbin.method, logbin.args)    
    mres <- match(outnames, names(res), 0L)
    fit[names(res)[mres]] <- res[mres]
    fit$family <- family
    good <- rep(TRUE, NROW(Y))
    w <- sqrt((fit$prior.weights[good] * family$mu.eta(fit$linear.predictors)[good]^2) / family$variance(fit$fitted.values)[good])
    z <- (fit$linear.predictors - (if(is.null(offset)) rep.int(0, NROW(Y)) else offset))[good] + 
      (fit$y - fit$fitted.values)[good] / family$mu.eta(fit$linear.predictors)[good]
    qr.tol <- min(1e-07, control$epsilon/1000)
    qr.val <- qr.default(fit$x * w, qr.tol, LAPACK = FALSE)
    fit$rank <- qr.val$rank
    qr.val$qr <- as.matrix(qr.val$qr)
    nvars <- ncol(fit$x)
    nr <- min(sum(good), nvars)
    if (nr < nvars) {
      Rmat <- diag(nvars)
      Rmat[1L:nr, 1L:nvars] <- qr.val$qr[1L:nr, 1L:nvars]
    }
    else Rmat <- qr.val$qr[1L:nvars, 1L:nvars]
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    fit$qr <- structure(append(qr.val[c("qr", "rank", "qraux", "pivot")], list(tol = qr.tol)), class = "qr")
    fit$effects <- qr.qty(qr.val, z * w)
    xnames <- dimnames(fit$x)[[2L]]
    xxnames <- xnames[fit$qr$pivot]
    colnames(fit$qr$qr) <- xxnames
    dimnames(Rmat) <- list(xxnames, xxnames)
    fit$R <- Rmat
    names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("", sum(good) - fit$rank))
    fit$weights <- rep.int(0, NROW(Y))
    fit$weights[good] <- w^2
    names(fit$weights) <- names(fit$y)
    if(model) fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    fit$terms <- mt
    fit$data <- data
    fit$offset <- offset
    fit$contrasts <- attr(fit$x, "contrasts")
    fit$xlevels <- .getXlevels(mt, mf)
    fit$xminmax <- .getXminmax(mt, mf)
  }
    
  fit$call <- call
  fit$formula <- formula
  fit$control <- control
  fit$method <- method
                    
  class(fit) <- c("logbin", "glm", "lm")
  fit
}