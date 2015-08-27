\name{conv.test}
\alias{conv.test}
\title{
Convergence Test Based on L2 Norm
}
\description{
Performs a test of convergence based on the L2 norm of
the change in the parameter estimates.
}
\usage{
conv.test(theta1, theta2, epsilon)
}
\arguments{
  \item{theta1}{
	vector of parameter estimates at previous
  step.
}
  \item{theta2}{
	vector of parameter estimates at current
  step.
}
  \item{epsilon}{
	positive convergence tolerance.
}
}
\details{
This is used as the convergence test in the
\code{\link{logbin}} fitting functions, because the EM
algorithm may converge slowly such that the test based on
the deviance used in \code{\link{glm.fit}} (see
\code{\link{glm.control}}) may report convergence at a
point away from the actual optimum.
}
\value{
A logical; \code{TRUE} if
\code{sqrt(sum((theta1-theta2)**2))/sqrt(sum(theta1**2)) <
epsilon}; \code{FALSE} otherwise.
}
\author{
Mark W. Donoghoe \email{mark.donoghoe@mq.edu.au}
}
\examples{
theta.old <- c(-4,-5,-6)
theta.new <- c(-4.05,-5,-6)

conv.test(theta.old, theta.new, 0.01)
conv.test(theta.old, theta.new, 0.005)
}
\keyword{ models }
\keyword{ misc }
