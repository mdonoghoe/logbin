logbin.design <- function(terms, data, type = c("cem", "em"), mats, rm = NULL) {
  type <- match.arg(type)
  X.orig <- model.matrix(terms, data)
  if (type == "em") {
    Amat.0 <- lapply(mats$Amat, function(x) x[1, , drop = FALSE])
    Amat.1 <- lapply(mats$Amat, function(x) x[-1, , drop = FALSE])
  } else {
    if (length(rm) != length(mats$Amat))
      stop("Length of reference vector is not equal to the number of terms")
    Amat.0 <- mapply(function(x, y, z) if (z > 1) x[1, -y, drop = FALSE] else 
                                              x[1, , drop = FALSE], mats$Amat, 
                     rm, mats$nref, SIMPLIFY = FALSE)
    Amat.1 <- mapply(function(x, y, z) if (z > 1) x[-1, -y, drop = FALSE] else
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