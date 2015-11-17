logbin.allref <- function(object, data = environment(object), type = c("cem","em"), mono, start = NULL) {
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
      stop("model frame and formula mismatch in logbin.allref()")
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
  
  nobs <- nrow(data)
  termlist <- attr(t, "term.labels")
  nvar <- length(termlist)
    
  npar <- sum(as.numeric(!isF))
  if(any(isF)) npar <- npar + sum(sapply(data[isF], function(x) nlevels(factor(x)) - 1)) 
  
  if (missing(mono)) mono <- rep(FALSE, nvar)
  if (is.null(mono)) mono <- rep(FALSE, nvar)
  monotonic <- rep(FALSE, nvar)
  names(monotonic) <- termlist
  monotonic[mono] <- TRUE
  names(monotonic) <- termlist
  
  allref <- list()
    
  if (!is.null(start)) {
    if (length(start) != npar)
      stop(gettextf("number of values in 'start' is %d should equal %d (number of parameters)",
           length(start), npar), domain = NA)
    start.orig <- start
    start.new.int <- start.orig[1]
    start.new.other <- start.orig[-1]
    this.start.o <- this.start.n <- 2
    delta.denom <- 1
  } else {
    start.new.int <- NULL
    start.new.other <- NULL
  }
  
  if (nvar == 0) return(list(allref = allref, terms = t, data = data))
  for (term in termlist) {
    allref[[term]] <- list()
    term2 <- gsub("`", "", term)
    if (!isF[term2]) {
      cont.min <- min(data[[term2]])
      cont.max <- max(data[[term2]])
      if (type == "em" || monotonic[term]) allref[[term]][[1]] <- 2
      else if (is.null(start)) allref[[term]] <- as.list(1:2)
      else allref[[term]] <- as.list(abs(3 * as.numeric(start.orig[this.start.o] > 0) - 1:2))
      if (!is.null(start)) {
        if (monotonic[term]) {
          start.new.int <- start.new.int + start.orig[this.start.o] * cont.max
          start.new.other[this.start.n-1] <- -start.orig[this.start.o]
        } else {
          cont.delta <- as.numeric(start.orig[this.start.o] <= 0)
          start.new.int <- start.new.int + start.orig[this.start.o] * (cont.delta * cont.min + (1 - cont.delta) * cont.max)
          if (type == "cem") start.new.other[this.start.n-1] <- cont.delta * start.orig[this.start.o] - (1 - cont.delta) * start.orig[this.start.o]
          else {
            start.new.other[this.start.n-1] <- -(1 - cont.delta) * start.orig[this.start.o]
            start.new.other <- c(start.new.other, cont.delta * start.orig[this.start.o])
          }
        }
        delta.denom <- delta.denom + cont.max - cont.min
        this.start.o <- this.start.o + 1
        this.start.n <- this.start.n + 1
      }
      attr(allref[[term]], "type") <- 1
    } else {
      lvls <- levels(factor(data[[term]]))
      nlvls <- nlevels(factor(data[[term]]))
      if (!is.null(start)) start.this <- start.orig[this.start.o:(this.start.o + nlvls - 2)]
      if (!monotonic[term]) {
        if (type == "em") allref[[term]][[1]] <- lvls[1]
        else if (is.null(start)) allref[[term]] <- as.list(rev(lvls))
        else allref[[term]] <- as.list(lvls[order(c(0, start.this), decreasing = TRUE)])
        if (!is.null(start)) {
          start.new.int <- start.new.int + max(c(0, start.this))
          start.new.other[(this.start.n - 1):(this.start.n + nlvls - 3)] <- (c(0, start.this) - max(c(0, start.this)))[-which.max(c(0, start.this))]
          if (type == "em") start.new.other <- append(start.new.other, 0, after = this.start.n - 3 + which.max(c(0, start.this)))
          this.start.o <- this.start.o + nlvls - 1
          this.start.n <- this.start.n + nlvls - as.numeric(type == "cem")
          delta.denom <- delta.denom + 1
        }
        attr(allref[[term]], "type") <- 2
      } else {
        allref[[term]][[1]] <- lvls
        if(!is.null(start)) {
          start.new.int <- start.new.int + start.this[length(start.this)]
          start.new.other[(this.start.n - 1):(this.start.n + nlvls - 3)] <- rev(diff(rev(c(0, start.this))))
          this.start.o <- this.start.o + nlvls - 1
          this.start.n <- this.start.n + nlvls - 1
        }
        attr(allref[[term]], "type") <- 3
      }
    }
  }  
  if(type == "em" && !is.null(start)) {
    start.new.int <- start.new.int / delta.denom
    start.new.other <- start.new.other + start.new.int
  }
  list(allref = allref, terms = t, data = data, monotonic = monotonic, start.new = c(start.new.int, start.new.other))
}