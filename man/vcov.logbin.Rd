\name{vcov.logbin}
\alias{vcov.logbin}
\title{
Calculate Variance-Covariance Matrix for a Fitted logbin Model Object
}
\description{
Returns the variance-covariance matrix of the main parameters of a fitted \code{logbin} model object.
}
\usage{
\method{vcov}{logbin}(object, ...)
}
\arguments{
  \item{object}{
	an object of class \code{"logbin"}, usually from a call to \code{\link{logbin}}
	or \code{\link{logbin.smooth}}.
}
  \item{\dots}{
    additional arguments for method functions.
}
}
\details{
An equivalent method to \code{\link{vcov}}, to use with \code{\link{logbin}} models.
}
\value{
A matrix of the estimated covariances between the parameter estimates in the linear or 
non-linear predictor of the model. This should have row and column names corresponding 
to the parameter names given by the \code{\link{coef}} method.
}
\author{
Mark W. Donoghoe \email{markdonoghoe@gmail.com}
}
\note{
If \code{object$boundary == TRUE}, the standard errors of the coefficients
are not valid, and a matrix of \code{NaN}s is returned.
}
\seealso{
\code{\link{summary.logbin}}, \code{\link{vcov.glm}}
}
\examples{
## For an example see example(logbin)
}
\keyword{ models }
\keyword{ regression }
