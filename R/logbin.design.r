logbin.design2 <- function(terms, data, type = c("cem", "em"), mats, rm = NULL) {
  type <- match.arg(type)
  X.orig <- model.matrix(terms, data)
  if (type == "em") {
    Amat.0 <- lapply(mats$Amat, function(x) x[1, , drop = FALSE])
    Amat.1 <- lapply(mats$Amat, function(x) x[-1, , drop = FALSE])
  } else {
    if (length(rm) != length(mats$Amat))
      stop("Length of reference vector is not equal to the number of terms")
    Amat.0 <- mapply(function(x, y, z) if (z > 1) x[1, -(ncol(x) - y + 1), drop = FALSE] else 
                                              x[1, , drop = FALSE], mats$Amat, 
                     rm, mats$nref, SIMPLIFY = FALSE)
    Amat.1 <- mapply(function(x, y, z) if (z > 1) x[-1, -(ncol(x) - y + 1), drop = FALSE] else
                                              x[-1, , drop = FALSE], mats$Amat, 
                     rm, mats$nref, SIMPLIFY = FALSE)
  }
  A0 <- do.call("cbind", Amat.0)
  A1 <- Matrix::.bdiag(Amat.1)
  A <- rbind(cbind(1, A0), cbind(0, A1))
  X.reparam <- cbind(1, X.orig %*% rbind(A0, A1))
  # Fix up possible numerical issues
  X.reparam[X.reparam <= .Machine$double.eps] <- 0
  list(X.orig = X.orig, X.reparam = X.reparam, A = A)
}

logbin.design <- function(terms, data, type = c("cem","em"), allref, mono, design.ref) {
  type <- match.arg(type)
  X.orig <- model.matrix(terms, data)
  terms.new <- terms
  data.new <- data
  termlabels <- attr(terms, "term.labels")
  design.type <- sapply(allref, attr, "type")
  ref.vector <- as.vector(design.ref, mode = "integer")
  for (i in seq_len(length(design.ref))) {
    varname <- gsub("`", "", termlabels[i])
    varref <- allref[[termlabels[i]]][[ref.vector[i]]]
    if (design.type[i] == 1) {
      cont.min <- min(data[[varname]], na.rm = TRUE)
      cont.max <- max(data[[varname]], na.rm = TRUE)
      if (varref == 1) data.new[[varname]] <- data[[varname]] - cont.min
      else data.new[[varname]] <- cont.max - data[[varname]]
      if (type == "em" && !mono[termlabels[i]]) {
        varname.rev <- paste0(varname, ".rev")
        if (varref == 1) data.new[[varname.rev]] <- cont.max - data[[varname]]
        else data.new[[varname.rev]] <- data[[varname]] - cont.min
        terms.new <- terms(update(terms.new, as.formula(paste(".~. +", varname.rev))))
      }
    } else if(design.type[i] == 2) {
      data.new[[varname]] <- relevel(factor(data[[varname]]), ref = varref)
      if(type == "cem")
        contrasts(data.new[[varname]]) <- contr.treatment(levels(data.new[[varname]]), base = 1)
      else
        contrasts(data.new[[varname]], nlevels(data.new[[varname]])) <- contr.treatment(levels(data.new[[varname]]), base = 1, contrasts = FALSE)
    } else if(design.type[i] == 3) {
      data.new[[varname]] <- factor(data[[varname]])
      contrasts(data.new[[varname]]) <- contr.isotonic.rev(levels(data.new[[varname]]), perm = varref)
    }
  }
  attr(data.new, "terms") <- terms.new
  X <- model.matrix(terms.new, data.new)  
  X
}