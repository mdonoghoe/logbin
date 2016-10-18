\name{logbin.smooth}
\alias{logbin.smooth}
\title{
Smooth Log-Binomial Regression
}
\description{
\code{logbin.smooth} fits log-link binomial
regression models using a stable CEM algorithm. It provides additional
flexibility over \code{\link{logbin}} by allowing for smooth 
semi-parametric terms.
}
\usage{
logbin.smooth(formula, mono = NULL, data, subset, na.action, offset, 
              control = list(...), model = TRUE, model.logbin = FALSE, 
              method = c("cem", "em"), accelerate = c("em", "squarem", "pem", "qn"),
              control.accelerate = list(), ...)
}
\arguments{
  \item{formula}{
	an object of class \code{"\link{formula}"}
  (or one that can be coerced into that class): a symbolic
  description of the model to be fitted. The details of
  model specification are given under "Details".  The model must contain an intercept
  and at least one semi-parametric term, included by using the
  \code{\link{B}} or \code{\link{Iso}} functions. Note that 2nd-order terms
  (such as interactions) or above are not currently supported
  (see \code{\link{logbin}}).
}
  \item{mono}{
	a vector indicating which terms in
  \code{formula} should be restricted to have a
  monotonically non-decreasing relationship with the
  outcome. May be specified as names or indices of the
  terms.
  
	\code{\link{Iso}()} terms are always monotonic.
}
  \item{data}{
	an optional data frame, list or environment
  (or object coercible by \code{\link{as.data.frame}} to a
  data frame) containing the variables in the model. If not
  found in \code{data}, the variables are taken from
  \code{\link{environment}(formula)}, typically the environment
  from which \code{logbin.smooth} is called.
}
  \item{subset}{
	an optional vector specifying a subset of
  observations to be used in the fitting process.
}
  \item{na.action}{
  a function which indicates what should happen when the data
  contain \code{NA}s. The default is set be the \code{na.action}
  setting of \code{\link{options}}, and is \code{\link{na.fail}}
  if that is unset. The `factory-fresh' default is \code{\link{na.omit}}.
  Another possible value is \code{NULL}, no action. Value
  \code{\link{na.exclude}} can be useful.
}
  \item{offset}{
	this can be used to specify an \emph{a
  priori} known component to be included in the linear
  predictor during fitting. This should be \code{NULL} or a
  \emph{non-positive} numeric vector of length equal to the number of cases.
  One or more \code{\link{offset}} terms can be included in
  the formula instead or as well, and if more than one is
  specified their sum is used. See
  \code{\link{model.offset}}.
}
  \item{control}{
	a list of parameters for controlling the
  fitting process, passed to
  \code{\link{logbin.control}}.
}
  \item{model}{
	a logical value indicating whether the \emph{model frame} 
	should be included as a component of the returned value.
}
  \item{model.logbin}{
	a logical value indicating whether the fitted \code{logbin}
	object should be included as a component of the returned value.
}
  \item{method}{
    a character string that determines which EM-type algorithm to use
  to find the MLE: \code{"cem"} for the combinatorial EM algorithm,
  which cycles through a sequence of constrained parameter spaces, or
  \code{"em"} for a single EM algorithm based on an overparameterised
  model.
  
  Unlike \code{\link{logbin}}, methods \code{"glm"} and \code{"ab"} are
  not available because they do not support the necessary monotonicity constraints.
}
  \item{accelerate}{
    a character string that determines the acceleration
  algorithm to be used, (partially) matching one of \code{"em"} (no acceleration -- the default), 
  \code{"squarem"}, \code{"pem"} or \code{"qn"}. See \code{\link[turboEM]{turboem}}
  for further details. Note that \code{"decme"} is not permitted.
}
  \item{control.accelerate}{
    a list of control parameters for the acceleration algorithm. See \code{\link[turboEM]{turboem}}
  for details of the parameters that apply to each algorithm. If not specified, the defaults are used.
}
  \item{\dots}{
	arguments to be used to form the default
  \code{control} argument if it is not supplied directly.
}
}
\details{
\code{logbin.smooth} performs the same fitting process as \code{\link{logbin}},
providing a stable maximum likelihood estimation procedure for log-link
binomial GLMs, with the added flexibility of allowing semi-parametric
\code{\link{B}} and \code{\link{Iso}} terms (note that \code{logbin.smooth} will stop with an 
error if no semi-parametric terms are specified in the right-hand side of the \code{formula}; 
\code{\link{logbin}} should be used instead).

The method partitions the parameter space associated with the semi-parametric part of the
model into a sequence of constrained parameter spaces, and defines a fully parametric 
\code{logbin} model for each. The model with the highest log-likelihood is the MLE for 
the semi-parametric model (see Donoghoe and Marschner, 2015).
}
\value{
An object of class \code{"logbin.smooth"}, which contains the same objects as class
\code{"logbin"} (the same as \code{"\link{glm}"}), as well as:
\item{model.logbin}{if \code{model.logbin} is \code{TRUE}; the \code{logbin} object
for the fully parametric model corresponding to the fitted model.}
\item{xminmax.smooth}{the minimum and maximum observed values for
each of the smooth terms in the model, to help define the 
covariate space.}
\item{full.formula}{the component from \code{\link{interpret.logbin.smooth}(formula)}
that contains the \code{formula} term with any additional arguments to the \code{\link{B}}
function removed.}
\item{knots}{a named list containing the knot vectors for each of the smooth
terms in the model.}
}
\references{
Donoghoe, M. W. and I. C. Marschner (2015). Flexible regression models for rate 
differences, risk differences and relative risks. \emph{International Journal of 
Biostatistics} 11(1): 91--108.

Marschner, I. C. (2014). Combinatorial EM algorithms.
\emph{Statistics and Computing} 24(6): 921--940.
}
\author{
Mark W. Donoghoe \email{markdonoghoe@gmail.com}
}
\seealso{
\code{\link{logbin}}
}
\examples{
## Simple example
x <- c(0.3, 0.2, 0.0, 0.1, 0.2, 0.1, 0.7, 0.2, 1.0, 0.9)
y <- c(5, 4, 6, 4, 7, 3, 6, 5, 9, 8)
system.time(m1 <- logbin.smooth(cbind(y, 10-y) ~ B(x, knot.range = 0:2), mono = 1, trace = 1))
## Compare with accelerated version
system.time(m1.acc <- update(m1, accelerate = "squarem"))
## Isotonic relationship
m2 <- logbin.smooth(cbind(y, 10-y) ~ Iso(x))
\donttest{
plot(m1)
plot(m2)
}
summary(predict(m1, type = "response"))
summary(predict(m2, type = "response"))
}
\keyword{ regression }
\keyword{ smooth }
