\name{logbin.allref}
\alias{logbin.allref}
\title{
Parameterisation for EM-type Algorithm
}
\description{
A workhorse function for \code{\link{logbin}},
\code{logbin.allref} takes the formula and data for a
log-link binomial GLM and produces a list of all
parameterisations needed for the associated CEM algorithm.
}
\usage{
logbin.allref(object, data = environment(object), type = c("cem", "em"), mono, start = NULL)
}
\arguments{
  \item{object}{
	a model \code{\link{formula}} or a
  \code{\link{terms}} object for the \code{\link{logbin}} model.
}
  \item{data}{
	a data frame created with
  \code{\link{model.frame}}. If another sort of object,
  \code{\link{model.frame}} is called first.
}
  \item{type}{
    a character string, indicating which type of parameterisation
  will be used. See ``Details" for more information.
}
  \item{mono}{
	a vector indicating which terms should be
  restricted to have a monotonically non-decreasing
  relationship with the outcome.
}
  \item{start}{
	starting values for the parameters in the linear predictor.
}
}
\details{
For a CEM algorithm (\code{type = "cem"}), the
parameter space is partitioned into a collection of
restricted parameter spaces (see Marschner, 2014).
\code{logbin.allref} finds the list of possible
parameterisations of each term in the model. 

If a term \code{x} has a \code{TRUE} value for
\code{\link{is.factor}(x)}, \code{\link{is.character}(x)}
or \code{\link{is.logical}(x)}, it is considered to be a
categorical covariate. This has a
parameterisation for each level of the factor.

Otherwise the covariate is considered to be continuous, in
which case it has two possible parameterisations, relating
to the minimum and maximum observed values.

If a covariate is restricted to be monotonic via the
\code{mono} argument, it has only one parameterisation.

\code{\link{logbin.cem}} considers all possible combinations of
the parameterisations of each covariate, and for each calls
\code{\link{logbin.design}} to create the appropriate
non-negative design matrix to be used in the EM algorithm.

For the overparameterised EM method (\code{type = "em"}),
only one parameter space is searched; this function
produces appropriate output that can be used with
\code{\link{logbin.design}}.
}
\value{
A list with components: 
\item{allref}{a named list, with
one component for each term in the model. Each component is
itself a list, whose components are each of the
parameterisations for that term. 

If \code{start} was specified, the first component for each
term will correspond to the parameterisation specified
by \code{start}.} 
\item{terms}{the \code{\link{terms}} component of \code{object}.}
\item{data}{the object passed into the \code{data}
argument, or the result of calling \code{\link{model.frame}} with \code{data}.}
\item{monotonic}{a named logical vector indicating which
components of \code{terms} are restricted to be
monotonically non-decreasing.}
\item{start.new}{a reparameterised version of \code{start}, corresponding
to the first parameterisation in \code{allref}. \code{NULL} if \code{start}
was not supplied.}
}
\references{
Marschner, I. C. (2014). Combinatorial EM algorithms.
\emph{Statistics and Computing} 24(6): 921--940.

Marschner, I. C. and A. C. Gillett (2012). Relative risk regression: reliable
and flexible methods for log-binomial models. \emph{Biostatistics} 13(1): 179--192.
}
\author{
Mark W. Donoghoe \email{mark.donoghoe@mq.edu.au}
}
\seealso{
\code{\link{logbin}}
}
\keyword{ misc }