logbin <- function (formula, mono = NULL, data, subset, na.action, start = NULL, offset,
            control = list(...), model = TRUE, method = c("cem","em","glm","glm2"),
            accelerate = c("em","squarem","pem","qn"), control.accelerate = list(list()), 
            warn = TRUE, ...) {
  call <- match.call()
  method <- match.arg(method)
  accelerate <- match.arg(accelerate)
  
  family <- binomial(link = log)
  if(missing(data)) data <- environment(formula)
  control <- do.call("logbin.control", control)
    
  outnames <- c("coefficients", "residuals", "fitted.values", "rank", "family",
                "linear.predictors", "deviance", "loglik", "aic", "aic.c",
                "null.deviance", "iter", "prior.weights", "weights",
                "df.residual", "df.null", "y", "x")
  if (model) outnames <- c(outnames, "model")
  outnames <- c(outnames, "converged", "boundary", "na.action", "call",
                "formula", "terms", "data", "offset", "control", "method", "xlevels",
                "xminmax", "np.coefficients", "nn.x")
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
  }
  else {
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
  
    logbin.method <- paste("logbin",method, sep = ".")
    res <- do.call(logbin.method, list(mt = mt, mf = mf, Y = Y, offset = offset, mono = mono,
                                       start = start, control = control, accelerate = accelerate, 
                                       control.accelerate = control.accelerate,
                                       warn = warn))    
    mres <- match(outnames, names(res), 0L)
    fit[names(res)[mres]] <- res[mres]
    fit$family <- family
    fit$weights <- rep(1, NROW(Y))
    if(model) fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    fit$terms <- mt
    fit$data <- data
    fit$offset <- offset
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