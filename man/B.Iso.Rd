\name{B.Iso}
\alias{B}
\alias{Iso}
\title{
Defining Smooths in logbin.smooth Formulae
}
\description{
Function used in the definition of smooth terms within \code{\link{logbin.smooth}}
model formulae. The function does not evaluate a smooth --- it exists purely to help
set up a model using smooths.
}
\usage{
B(..., knots = NULL, knot.range = 0:5)

Iso(...)
}
\arguments{
  \item{\dots}{
	variable that this smooth is a function of. Note that unlike \code{\link[mgcv]{gam}},
	smooths that are functions of more than one variable are not supported.
}
  \item{knots}{
	\emph{unique} positions of \emph{interior} knots of a B-spline basis. Boundary
	knots are created automatically.
}
  \item{knot.range}{
	if \code{knots} is not specified, a vector containing a series of non-negative
	integers denoting the number of \emph{interior} knots for which the model will be fit.
	These are placed at evenly-spaced quantiles of the observed covariate values.
	
	At least one of \code{knots} or \code{knot.range} must be non-missing.
}
}
\details{
The function does not evaluate the variable arguments; the
output from this function is used when producing the model
matrix, at which point the actual basis functions are constructed.

\code{B} is used to specify an order-3 B-spline basis (which can be
restricted to be monotonically non-decreasing via the
\code{mono} argument in \code{\link{logbin.smooth}}). If
\code{length(knot.range) > 1}, models with each of the
specified number of interior knots will be fit, and the model
with the best (smallest) \code{aic.c} will be returned.

\code{Iso} is used to specify an isotonic basis, designed
such that the resulting function has non-negative
increments at each observed covariate value. When
\code{Iso} is used, the resulting function will always be
monotonically non-decreasing, regardless of the value of
\code{mono}.
}
\value{
An object of class \code{"B.smooth"} (for \code{B}) or
\code{"Iso.smooth"} (for \code{Iso}), which is a list with
the following elements: \item{term}{name of the term
provided in the \code{...} argument.}
\item{termlabel}{label for the term in the model; e.g. for
term \code{"x"} it will be \code{"B(x)"} or
\code{"Iso(x)"}.} \item{knots}{vector of interior knots (if
specified). \code{NA} for \code{Iso}.}
\item{knot.range}{vector of number of interior knots.
\code{NA} for \code{Iso}.}
}
\author{
Mark W. Donoghoe \email{mark.donoghoe@mq.edu.au}
}
\seealso{
\code{\link{logbin.smooth}}

\code{\link[mgcv]{s}} performs a similar function in the \code{\link[mgcv]{mgcv}} package.
}
\examples{
## See example(logbin.smooth) for an example of specifying smooths in model
##  formulae.
}
\keyword{ smooth }
