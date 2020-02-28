# Produces a convergence test function with default values as supplied
conv.test <- function(Amat, epsilon) {
  conv.user <- function(old, new, Amat, epsilon) {
    old.r <- as.matrix(Amat %*% old)
    new.r <- as.matrix(Amat %*% new)
    diff.r <- sqrt(crossprod(old.r - new.r) / crossprod(old.r))
    (diff.r < epsilon)
  }
  formals(conv.user)$Amat <- Amat
  formals(conv.user)$epsilon <- epsilon
  conv.user
}