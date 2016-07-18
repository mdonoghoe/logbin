\name{summary.logbin}
\alias{summary.logbin}
\alias{print.summary.logbin}
\title{
Summarising logbin Model Fits
}
\description{
These functions are all \code{\link{methods}} for class \code{logbin}
or \code{summary.logbin} objects.
}
\usage{
\method{summary}{logbin}(object, correlation = FALSE, ...)

\method{print}{summary.logbin}(x, digits = max(3L, getOption("digits") - 3L), 
      signif.stars = getOption("show.signif.stars"), ...)
}
\arguments{
  \item{object}{
	an object of class \code{"logbin"}, usually from a call to \code{\link{logbin}}
	or \code{\link{logbin.smooth}}.
}
  \item{x}{
	an object of class \code{"summary.logbin"}, usually from a call to
	\code{summary.logbin}.
}
  \item{correlation}{
	logical; if \code{TRUE}, the correlation matrix of the estimated parameters is
	returned and printed.
}
  \item{digits}{
	the number of significant digits to use when printing.
}
  \item{signif.stars}{
	logical; if \code{TRUE}, `significance stars' are printed for each coefficient.
}
  \item{\dots}{
	further arguments passed to or from other methods.
}
}
\details{
These perform the same function as \code{\link{summary.glm}} and \code{\link{print.summary.glm}},
producing similar results for \code{logbin} models. \code{print.summary.logbin} additionally prints 
the small-sample corrected AIC (\code{aic.c}), and the number of EM iterations for the parameterisation 
corresponding to the MLE.

The dispersion used in calculating standard errors is fixed as \code{1}.
}
\value{
\code{summary.logbin} returns an object of class \code{"summary.logbin"}, a list with components
\item{call}{the component from \code{object}.}
\item{family}{the component from \code{object}.}
\item{deviance}{the component from \code{object}.}
\item{aic}{the component from \code{object}.}
\item{aic.c}{the component from \code{object}.}
\item{df.residual}{the component from \code{object}.}
\item{null.deviance}{the component from \code{object}.}
\item{df.null}{the component from \code{object}.}
\item{iter}{the component from \code{object}.}
\item{deviance.resid}{the deviance residuals: see \code{\link{residuals.glm}}.}
\item{coefficients}{the matrix of coefficients, standard errors, z-values and p-values.}
\item{aliased}{included for compatibility --- always \code{FALSE}.}
\item{dispersion}{the inferred/estimated dispersion.}
\item{df}{included for compatibility --- a 3-vector of the number of coefficients,
the number of residual degrees of freedom, and the number of coefficients (again).}
\item{cov.unscaled}{the unscaled (\code{dispersion = 1}) estimated covariance
matrix of the estimated coefficients. \code{NaN} if \code{object$boundary == TRUE}.}
\item{cov.scaled}{ditto, scaled by \code{dispersion}.}
\item{correlation}{if \code{correlation} is \code{TRUE}, the estimated correlations
of the estimated coefficients. \code{NaN} if \code{object$boundary == TRUE}.}
}
\author{
Mark W. Donoghoe \email{markdonoghoe@gmail.com}
}
\note{
If \code{object$boundary == TRUE}, the standard errors of the coefficients
are not valid, and a matrix of \code{NaN}s is returned by \code{\link{vcov.logbin}}. If
the MLE is not on the boundary but the model contains parameters with monotonicity 
constraints, the standard errors do not take this into account and should be used
with caution.
}
\seealso{
\code{\link{logbin}}, \code{\link{summary.glm}}
}
\examples{
## For examples see example(logbin)
}
\keyword{ models }
\keyword{ regression }
