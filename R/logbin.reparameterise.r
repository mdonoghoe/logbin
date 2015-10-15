logbin.reparameterise <- function(np.coefs, terms, data, type = c("cem","em"), allref, mono, design.ref)
{
  type <- match.arg(type)
  termlabels <- attr(terms, "term.labels")
  design.type <- sapply(allref, attr, "type")
  ref.vector <- as.vector(design.ref, mode = "integer")
    
  npar.o <- npar.n <- 1L + sum(design.type == 1)
  if(sum(design.type == 2) > 0) {
    npar.o <- npar.o + sum(sapply(sapply(termlabels[design.type == 2], function(x) gsub("`","",x)), 
                                  function(x) nlevels(factor(data[[x]])) - as.numeric(type == "cem")))
    npar.n <- npar.n + sum(sapply(sapply(termlabels[design.type == 2], function(x) gsub("`","",x)), 
                                  function(x) nlevels(factor(data[[x]])) - 1))
  }
  if(sum(design.type == 3) > 0) {
    npar.o <- npar.o + sum(sapply(sapply(termlabels[design.type == 3], function(x) gsub("`","",x)), 
                                  function(x) nlevels(factor(data[[x]])) - 1))
    npar.n <- npar.n + sum(sapply(sapply(termlabels[design.type == 3], function(x) gsub("`","",x)), 
                                  function(x) nlevels(factor(data[[x]])) - 1))
  }
    
  coefs.int <- coefs.int.reparam <- coefs.boundary <- np.coefs[1]
  coefs.model <- np.coefs[-1L]
  coefs.model.reparam <- rep(0, npar.n - 1)
  coef.count.o <- coef.count.n <- coef.cont.count <- 0L
  
  coef.names <- "(Intercept)"
    
  for (i in seq_len(length(design.ref))) {
    varname <- gsub("`", "", termlabels[i])
    varref <- allref[[termlabels[i]]][[ref.vector[i]]]
    if (design.type[i] == 1) {
      thiscoef <- coefs.model[(coef.count.o + 1L)]
      cont.min <- min(data[[varname]])
      cont.max <- max(data[[varname]])
      if (type == "cem" || mono[termlabels[i]]) {
        if (mono[termlabels[i]]) coefs.boundary <- append(coefs.boundary, thiscoef)
        if (varref == 1) {
            coefs.int.reparam <- coefs.int.reparam - thiscoef * cont.min
            coefs.model.reparam[(coef.count.n + 1L)] <- thiscoef
        } else {
            coefs.int.reparam <- coefs.int.reparam + thiscoef * cont.max
            coefs.model.reparam[(coef.count.n + 1L)] <- -thiscoef
        }
      } else {
        thiscoef.rev <- coefs.model[(npar.o + coef.cont.count)]
        coefs.int.reparam <- coefs.int.reparam + thiscoef * cont.max - thiscoef.rev * cont.min
        coefs.model.reparam[(coef.count.n + 1L)] <- thiscoef.rev - thiscoef
      }
      coef.names <- append(coef.names, varname)
      coef.count.o <- coef.count.o + 1L
      coef.count.n <- coef.count.n + 1L
      coef.cont.count <- coef.cont.count + as.numeric(!mono[termlabels[i]])
    } else if (design.type[i] == 2) {
      lev <- levels(factor(data[[varname]]))
      nlev <- nlevels(factor(data[[varname]]))
      ref.orig <- varref
      ref.new <- lev[1L]
      if (type == "cem") {
        if (ref.orig != ref.new) {
          thiscoef <- coefs.model[(coef.count.o + 1L):(coef.count.o + nlev - 1L)]
          thiscoef.new <- append(thiscoef, 0, after = which(lev == ref.orig) - 1L)[-1L] - thiscoef[1L]
          coefs.int.reparam <- coefs.int.reparam + thiscoef[1L]
          coefs.model.reparam[(coef.count.n + 1L):(coef.count.n + nlev - 1L)] <- thiscoef.new
        }
      } else {
        thiscoef <- coefs.model[(coef.count.o + 1L):(coef.count.o + nlev)]
        thiscoef.new <- thiscoef[-1L] - thiscoef[1L]
        coefs.int.reparam <- coefs.int.reparam - max(c(0,thiscoef.new))
        coefs.boundary <- coefs.boundary + max(thiscoef)
        coefs.model.reparam[(coef.count.n + 1L):(coef.count.n + nlev - 1L)] <- thiscoef.new
      }
      coef.names <- append(coef.names, paste0(varname, lev[-1L]))
      coef.count.o <- coef.count.o + nlev - as.numeric(type == "cem")
      coef.count.n <- coef.count.n + nlev - 1
    } else if (design.type[i] == 3) {
      lev <- levels(factor(data[[varname]]))
      nlev <- nlevels(factor(data[[varname]]))
      ref.orig <- varref
      thiscoef <- coefs.model[(coef.count.o + 1L):(coef.count.o + nlev - 1L)]
      thiscoef.sum <- append(rev(cumsum(rev(thiscoef))), 0)
      thiscoef.ord <- thiscoef.sum[match(lev, ref.orig)]
      thiscoef.new <- thiscoef.ord[-1L] - thiscoef.ord[1L]
      coefs.int.reparam <- coefs.int.reparam + thiscoef.ord[1L]
      coefs.boundary <- append(coefs.boundary, thiscoef)
      coefs.model.reparam[(coef.count.n + 1L):(coef.count.n + nlev - 1L)] <- thiscoef.new
      coef.names <- append(coef.names, paste0(varname, lev[-1L]))
      coef.count.o <- coef.count.o + nlev - 1L
      coef.count.n <- coef.count.n + nlev - 1L
    }
  }
  coefs <- c(coefs.int.reparam, coefs.model.reparam)
  names(coefs) <- coef.names
  
  design <- model.matrix(terms, data)
  
  list(coefs = coefs, design = design, coefs.boundary = coefs.boundary)
}