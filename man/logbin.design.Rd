\name{logbin.design}
\alias{logbin.design}
\title{
Construct Design Matrix for logbin Model
}
\description{
\code{logbin.design} constructs the design matrix for a
\code{\link{logbin}} model, given a particular parameterisation.

This is a workhorse function --- it would not normally be called directly.
}
\usage{
logbin.design(terms, data, type = c("cem", "em"), allref, mono, design.ref)
}
\arguments{
  \item{terms}{
	\code{terms} component of object returned
  from a call to \code{\link{logbin.allref}} for the
  desired model.
}
  \item{data}{
	\code{data} component of object returned from
  a call to \code{\link{logbin.allref}} for the desired
  model.
}
  \item{type}{
    a character string, indicating which type of parameterisation
  will be used. See ``Details" for more information.
}
  \item{allref}{
	\code{allref} component of object returned
  from a call to \code{\link{logbin.allref}} for the
  desired model.
}
  \item{mono}{
    a logical vector, with names corresponding to the \code{"term.labels"}
  attribute of \code{terms}, indicating which terms should be
  restricted to have a monotonically non-decreasing
  relationship with the outcome
}
  \item{design.ref}{
	vector of indices for a particular
  parameterisation in \code{allref}. That is, each element
  corresponds to a term \code{x} in the model, and the
  value of the element indicates which item in the list
  \code{allref[[x]]} is the reference level in this
  parameterisation.
}
}
\details{
For a CEM algorithm (\code{type = "cem"}),
we must consider the
Cartesian product of all possible parameterisations. The
list of these for each term in the model is constructed by
a call to \code{\link{logbin.allref}}, and a list of all
possible combinations created using
\code{\link{expand.grid}}.

For a particular combination of reference levels,
\code{logbin.design} constructs the associated design
matrix by transforming \code{data}. Continuous covariates
are transformed such that their minimum or maximum observed
value corresponds to a transformed value of zero;
categorical covariates are transformed by using either
\code{\link{contr.treatment}} with a specified reference
level or \code{\link{contr.isotonic.rev}} so that the levels
are increasing in the specified order.

For the EM algorithm (\code{type = "em"}), the model
is overparameterised by introducing an additional parameter
for each non-monotonic categorical or continuous covariate.
\code{logbin.design} constructs the associated design matrix
by transforming \code{data} and including additional covariates
as appropriate.
}
\value{
A strictly non-negative design matrix to be passed to the
relevant fitting function.
}
\author{
Mark W. Donoghoe \email{mark.donoghoe@mq.edu.au}
}
\seealso{
\code{\link{logbin.allref}}, \code{\link{model.matrix}},
\code{\link{contr.treatment}}, \code{\link{contr.isotonic.rev}}
}
\keyword{ models }
