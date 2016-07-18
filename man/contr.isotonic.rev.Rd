\name{contr.isotonic.rev}
\alias{contr.isotonic.rev}
\title{
Contrast Matrix for Reversed Isotonic Covariate
}
\description{
Return something similar to a contrast matrix for a categorical
covariate that we wish to be monotonically non-decreasing
in a specified order.
}
\usage{
contr.isotonic.rev(n, perm, contrasts = TRUE, sparse = FALSE)
}
\arguments{
  \item{n}{
	a vector of levels for a factor, or the number
  of levels.
}
  \item{perm}{
	a permutation of the levels of \code{n} (or
  of the numbers \code{1:n}), which define the order in
  which the coefficients must be monotonically
  non-decreasing.
}
  \item{contrasts}{
	a logical indicating whether constrasts
  should be computed.
}
  \item{sparse}{
	included for compatibility reasons. Has no
  effect.
}
}
\details{
This function is used in creating the design matrix
for categorical covariates with a specified order under a
particular parameterisation. This is required
if a categorical covariate is defined as monotonic.

In the order specified by \code{perm}, the coefficient
associated with each level is the sum of increments between
the following levels. That is, if there are a total of \eqn{k}
levels, the first level is defined as \eqn{d_2 + d_3 + d_4 + \cdots + d_k}, 
the second as \eqn{d_3 + d_4 + \cdots + d_k}, 
the third as \eqn{d_4 + \cdots + d_k}, and so on.  In fitting the model, 
these increments are constrained to be non-positive.

Note that these are not `contrasts' as defined in the
theory for linear models, rather this is used to define the
\code{contrasts} attribute of each variable so that
\code{\link{model.matrix}} produces the desired design
matrix.
}
\value{
A matrix with \code{n} rows and \eqn{k} columns, with
\eqn{k=n-1} if \code{contrasts} is \code{TRUE} and
\eqn{k=n} if \code{contrasts} is \code{FALSE}.
}
\author{
Mark W. Donoghoe \email{markdonoghoe@gmail.com}
}
\seealso{
\code{\link{model.matrix}}, which uses
\code{contr.isotonic.rev} to create the design matrix.

\code{\link{contr.treatment}}, \code{\link{contrasts}} for
their usual use in regression models.
}
\examples{
contr.isotonic.rev(4,1:4)
contr.isotonic.rev(4,c(1,3,2,4))

# Show how contr.isotonic.rev applies within model.matrix
x <- factor(round(runif(20,0,2)))
mf <- model.frame(~x)
contrasts(x) <- contr.isotonic.rev(levels(x), levels(x))
model.matrix(mf)
}
\keyword{ design }
