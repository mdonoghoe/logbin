\name{logbin.smooth.design}
\alias{logbin.smooth.design}
\title{
Construct Design Matrix for logbin.smooth Model
}
\description{
\code{logbin.smooth.design} constructs the design matrix and associated formula
for an \code{\link{logbin.smooth}} model, given a particular parameterisation,
to be passed into \code{\link{logbin}} for fitting.

This is a workhorse function --- it would not normally be called directly.
}
\usage{
logbin.smooth.design(interpret, type = c("cem", "em"), allref, design.knots, design.param)
}
\arguments{
  \item{interpret}{
	the object returned by running \code{\link{interpret.logbin.smooth}}
	for the desired model, containing details of the smooth components.
}
  \item{type}{
    a character string, indicating which type of parameterisation
  will be used. See ``Details" for more information.
}
  \item{allref}{
	the object returned by running \code{\link{logbin.smooth.allref}}
	for the desired model.
}
  \item{design.knots}{
	a vector containing the number of internal knots for each smooth
	term (\code{NA} for \code{\link{Iso}} terms).
}
  \item{design.param}{
	a vector of indices for a particular parameterisation in 
	\code{allref$allref}. Each element corresponds to a smooth
	term in the model, and the value indicates which item in
	the associated list is the reference level for this
	parameterisation.
}
}
\details{
For a particular combination of reference levels, \code{logbin.smooth.design}
constructs the associated design matrix and formula. Specifically, for \code{\link{Iso}} smooth
components, it creates the matrix of indicator covariates for increments between
levels. For \code{\link{B}} smooth components, it creates the basis functions using 
\code{\link[splines]{splineDesign}} and --- if \code{type = "cem"} --- removing the 
column associated with the reference level (see Donoghoe and Marschner, 2015).

The \code{formula} component is altered to include the terms in the design matrix,
and \code{allref$monotonic} is altered such that all of the smooth coefficients
are restricted to be non-negative, as required.
}
\value{
A list with components:
\item{formula}{an updated version of \code{interpret$full.formula}, with
smooth terms removed and replaced by the names of their associated
basis components.}
\item{data}{an updated version of \code{interpret$data}, with columns
for the basis functions of the smooth terms added.}
\item{monotonic}{an updated version of \code{allref$monotonic}, such that
the coefficients associated with the smooth terms for this parameterisation
are constrained to be non-negative.}
\item{knots}{a list, with one component for each smooth term, containing
the knot vector for that term (\code{NA} for \code{\link{Iso}} terms).}
}
\references{
Donoghoe, M. W. and I. C. Marschner (2015). Flexible regression models for rate 
differences, risk differences and relative risks. \emph{International Journal of 
Biostatistics} 11(1): 91--108.
}
\author{
Mark W. Donoghoe \email{mark.donoghoe@mq.edu.au}
}
\seealso{
\code{\link{logbin.smooth}}, \code{\link{interpret.logbin.smooth}},
\code{\link{logbin.smooth.allref}}, \code{\link{logbin.design}}
}
\keyword{ models }
\keyword{ smooth }
