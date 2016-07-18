\name{logbin.control}
\alias{logbin.control}
\title{
Auxiliary for Controlling logbin Fitting
}
\description{
Auxiliary function for \code{\link{logbin}} fitting.
Typically only used internally by \code{\link{nplbin}}, 
but may be used to construct a
\code{control} argument to that function.
}
\usage{
logbin.control(bound.tol = 1e-06, epsilon = 1e-08, maxit = 10000, trace = 0)
}
\arguments{
  \item{bound.tol}{
	positive tolerance specifying the
  interior of the parameter space. If the fitted model is
  more than \code{bound.tol} away from the boundary of the
  parameter space then it is assumed to be in the interior.
  This can allow the computational method to terminate
  early if an interior maximum is found.  No early
  termination is attempted if \code{bound.tol = Inf}.
}
  \item{epsilon}{
	positive convergence tolerance
  \eqn{\epsilon}; the estimates are considered to have 
  converged when \eqn{\sqrt{
  \sum (\theta_{old} - \theta_{new})^2} / \sqrt {\sum
  \theta_{old}^2} < \epsilon}, where \eqn{\theta} is the
  vector of parameter estimates. See \code{\link{conv.test}}.
  
  This should be smaller than \code{bound.tol}.
}
  \item{maxit}{
	integer giving the maximum number of EM
  algorithm iterations for a given parameterisation.
}
  \item{trace}{
	number indicating level of output that
  should be produced. >= 1 gives output for each
  parameterisation, >= 2 gives output at each iteration.
}
}
\details{
This is used similarly to \code{\link{glm.control}}. The
\code{control} argument of \code{\link{logbin}} is by
default passed to the \code{control} argument of
\code{\link{nplbin}}.

When \code{trace} is greater than zero, calls to
\code{\link{cat}} produce the output. Hence,
\code{\link{options}(digits = *)} can be used to increase
the precision.
}
\value{
A list with components named as the arguments.
}
\author{
Mark W. Donoghoe \email{markdonoghoe@gmail.com}
}
\seealso{
\code{\link{glm.control}}, the equivalent function for
\code{\link{glm}} fitting.

\code{\link{nplbin}}, the
function used to fit \code{\link{logbin}} models.
}
\examples{
## Variation on example(glm.control) :

evts <- c(18,17,15,20,10,20,25,13,12)
obs <- rep(30,9)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
oo <- options(digits = 12)
logbin.D93X <- logbin(cbind(evts,obs-evts) ~ outcome + treatment, trace = 2, epsilon = 1e-2)
options(oo)
coef(logbin.D93X)
}
\keyword{ optimize }
\keyword{ models }
