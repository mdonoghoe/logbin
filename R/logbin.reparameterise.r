logbin.reparam <- function(object, data = environment(object), 
                           type = c("cem","em"), mono) {
  type <- match.arg(type)
  t <- if (missing(data))
    terms(object)
  else terms(object, data = data)
  if (is.null(attr(data, "terms")))
    data <- model.frame(object, data)
  else {
    reorder = match(sapply(attr(t, "variables"), deparse,
                           width.cutoff = 500)[-1L], names(data))
    if (any(is.na(reorder)))
      stop("model frame and formula mismatch in logbin.reparam()")
    if (!identical(reorder, seq_len(ncol(data))))
      data <- data[, reorder, drop = FALSE]
  }
  int <- attr(t, "response")
  
  namD <- names(data)
  for (i in namD) 
    if (is.character(data[[i]]))
      data[[i]] <- factor(data[[i]])
  isF <- vapply(data, function(x) is.factor(x) || is.logical(x), NA)
  isF[int] <- FALSE
  
  termlist <- attr(t, "term.labels")
  nvar <- length(termlist)
  
  if (missing(mono)) mono <- rep(FALSE, nvar)
  if (is.null(mono)) mono <- rep(FALSE, nvar)
  monotonic <- rep(FALSE, nvar)
  names(monotonic) <- termlist
  monotonic[mono] <- TRUE
  names(monotonic) <- termlist
  
  Vmat <- Amat <- nref <- npar <- list()
  
  for (term in termlist) {
    term2 <- gsub("`", "", term)
    if (!isF[term2]) {
      Vmat[[term]] <- matrix(range(data[[term2]]), ncol = 1)
      if (monotonic[term]) {
        rmult <- matrix(1:0, ncol = 1)
        nref[[term]] <- 1L
      } else {
        rmult <- diag(nrow(Vmat[[term]]))
        nref[[term]] <- nrow(Vmat[[term]])
      }
      npar[[term]] <- 1L
    } else {
      lvls <- levels(factor(data[[term]]))
      nlvls <- nlevels(factor(data[[term]]))
      Vmat[[term]] <- rbind(0, diag(nlvls-1))
      rmult <- diag(nlvls)
      nref[[term]] <- nlvls
      if (monotonic[term]) {
        rmult[upper.tri(rmult)] <- 1
        rmult <- rmult[,-ncol(rmult)]
        nref[[term]] <- 1
      }
      npar[[term]] <- nlvls - 1L
    }
    colnames(Vmat[[term]]) <- paste0(term2, ".s", seq(ncol(Vmat[[term]])))
    Amat[[term]] <- solve(cbind(1, Vmat[[term]])) %*% rmult
  }
  list(Vmat = Vmat, Amat = Amat, nref = nref, npar = npar,
       monotonic = monotonic, termlist = termlist, isF = isF)
}

# Go from transformed parameters (all non-positive), back to the original scale
logbin.reduce <- function(np.coefs, Amat) {
  as.matrix(Amat %*% np.coefs)
}

# Used when np.coefs is a matrix
coefhist.reduce <- function(np.coefhist, Amat, cnames) {
  coefhist <- as.matrix(t(logbin.reduce(t(np.coefhist), Amat)))
  colnames(coefhist) <- cnames
  rownames(coefhist) <- rownames(np.coefhist)
  coefhist
}

logbin.expand <- function(coefs, reparam, type = c("cem", "em")) {
  type <- match.arg(type)
  # Check that we have the right number of parameters
  nvar <- length(reparam$npar)
  nvar.nonmono <- sum(as.numeric(!reparam$monotonic))
  npar <- 1L + sum(unlist(reparam$npar))
  if (length(coefs) != npar)
    stop(gettextf("number of values in 'coefs' is %d should equal %d (number of parameters)",
                  length(coefs), npar), domain = NA)
  # Split coefs into covariates
  coefids <- rep(seq_along(reparam$termlist), unlist(reparam$npar))
  coefs.int <- coefs[1L]
  coefs.cov <- split(coefs[-1L], coefids)
  names(coefs.cov) <- reparam$termlist
  vstar.id <- astar.0 <- astar.1 <- list()
  for (term in reparam$termlist) {
    term2 <- gsub("`", "", term)
    coefs.this <- coefs.cov[[term]]
    if (reparam$monotonic[term]) {
      vstar.id[[term]] <- 1L
      astar.1[[term]] <- solve(reparam$Amat[[term]][-1L, , drop = FALSE])
      astar.0[[term]] <- -reparam$Amat[[term]][1L, , drop = FALSE] %*% astar.1[[term]]
    } else {
      # Identify the vertex that gives the largest fitted value
      vstar.id[[term]] <- which.max(c(0, coefs.this))
      vstar <- reparam$Vmat[[term]][vstar.id[[term]], , drop = FALSE]
      # Get the components of the expansion matrix
      astar.0[[term]] <- vstar
      astar.1[[term]] <- reparam$Vmat[[term]] - matrix(vstar, nrow = nrow(reparam$Vmat[[term]]),
                                                       ncol = ncol(vstar), byrow = TRUE)
      if (type == "cem") astar.1[[term]] <- astar.1[[term]][-vstar.id[[term]], , drop = FALSE]
    }
  }
  astar <- rbind(cbind(1, do.call("cbind", astar.0)),
                 cbind(0, Matrix::.bdiag(astar.1)))
  coefs.exp <- as.vector(astar %*% coefs)
  # Indicator for which 0 coefficients would correspond to boundary values
  # (the intercept and any constrained by monotonicity)
  npar.exp <- sapply(astar.1, nrow)
  which.boundary.split <- split(rep(FALSE, length(coefs.exp) - 1),
                                rep(seq_along(reparam$termlist), npar.exp))
  which.boundary.split[reparam$monotonic] <- lapply(which.boundary.split[reparam$monotonic], "!")
  which.boundary <- c(TRUE, unlist(which.boundary.split))
  # For type = "em", add something so that coefficients are not exactly 0
  if (type == "em" & nvar.nonmono > 0) {
    # Take the intercept as close to 0 as we can
    delta <- coefs.exp[1L] / (nvar.nonmono + 1)
    # Split coefs.exp into covariates
    coefids.exp <- rep(seq_along(reparam$termlist), unlist(reparam$npar) + as.numeric(!reparam$monotonic))
    coefs.int.exp <- coefs.exp[1L]
    coefs.cov.exp <- split(coefs.exp[-1L], coefids.exp)
    names(coefs.cov.exp) <- reparam$termlist
    # Only update the non-monotonic covariates
    coefs.int.exp <- delta
    coefs.cov.exp[!reparam$monotonic] <- lapply(coefs.cov.exp[!reparam$monotonic],
                                                function(x) x + delta)
    coefs.exp <- c(coefs.int.exp, unlist(coefs.cov.exp))
  }
  list(coefs.exp = coefs.exp, vstar.id = vstar.id, which.boundary = which.boundary)
}