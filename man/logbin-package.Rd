\name{logbin-package}
\alias{logbin-package}
\docType{package}
\title{
Relative Risk Regression Using the Log-Binomial Model
}
\description{
Methods for fitting log-link GLMs and GAMs to binomial data, using EM-type
algorithms with more stable convergence properties than standard methods.
}
\details{
\tabular{ll}{
Package: \tab logbin\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2015-10-19\cr
License: \tab GPL (>= 2)\cr
}
This package provides methods to fit generalised linear models (GLMs) and generalised additive
models (GAMs) with log link functions to binomial data.
It has two primary functions: \code{\link{logbin}} and \code{\link{logbin.smooth}}, 
together with various supporting functions.

It is useful in two main situations. The first is when a standard GLM routine, such as
\code{\link{glm}}, fails to converge with such a model. The second is when a flexible
semi-parametric component is desired in these models. One of the main purposes of this
package is to provide parametric and semi-parametric adjustment of relative risks.

The computational method, determined by the \code{method} argument, is either a 
combinatorial EM (CEM) algorithm (Marschner, 2014) or an overparameterised EM approach, 
which accommodate the parameter constraints and are more stable than iteratively
reweighted least squares. 

In a CEM algorithm, a collection of restricted parameter spaces is defined
which covers the full parameter space, and an EM algorithm is applied within each
restricted parameter space in order to find a collection of restricted maxima of
the log-likelihood function, from which can be obtained the global maximum over the
full parameter space.

In the overparameterised EM approach, additional parameters are added to the model,
and an EM algorithm finds the MLE of this overparameterised model by imposing
constraints on each individual parameter. This requires a single application of
the EM algorithm.

In each case, the EM algorithm may be accelerated by using the capabilities of the
\code{\link[turboEM]{turboem}} package.

The methodology implemented in this package is presented in Marschner and Gillett (2012) 
and Donoghoe and Marschner (2015).
}
\author{
Mark W. Donoghoe \email{mark.donoghoe@mq.edu.au}

Maintainer: Mark W. Donoghoe \email{mark.donoghoe@mq.edu.au}
}
\references{
Donoghoe, M. W. and I. C. Marschner (2015). Flexible regression models for rate 
differences, risk differences and relative risks. \emph{International Journal of 
Biostatistics} 11(1): 91--108.

Marschner, I. C. (2014). Combinatorial EM algorithms.
\emph{Statistics and Computing} 24(6): 921--940.

Marschner, I. C. and A. C. Gillett (2012). Relative risk regression: reliable
and flexible methods for log-binomial models. \emph{Biostatistics} 13(1): 179--192.
}
\keyword{ package }
\keyword{ regression }
\seealso{
\code{\link{glm}}, \code{\link[turboEM]{turboem}}
}
\examples{
## For examples, see example(logbin) and example(logbin.smooth)
}
