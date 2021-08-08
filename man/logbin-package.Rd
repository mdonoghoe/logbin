\name{logbin-package}
\alias{logbin-package}
\docType{package}
\title{
Relative Risk Regression Using the Log-Binomial Model
}
\description{
Methods for fitting log-link GLMs and GAMs to binomial data, including EM-type
algorithms with more stable convergence properties than standard methods.
}
\details{
\tabular{ll}{
Package: \tab logbin\cr
Type: \tab Package\cr
Version: \tab 2.0.5\cr
License: \tab GPL (>= 2)\cr
}
This package provides methods to fit generalised linear models (GLMs) and generalised additive
models (GAMs) with log link functions to binomial data, which can be used to
estimate adjusted relative risks. It has two primary functions: \code{\link{logbin}} 
and \code{\link{logbin.smooth}}, together with various supporting functions.

Standard GLM routines such as base R's \code{\link{glm}} typically use a modified Fisher
scoring algorithm, but this can experience numerical problems and fail to converge to the
maximum likelihood estimate (MLE). The \pkg{\link[=glm2-package]{glm2}}
package improves on this but can still have difficulties, particularly when the MLE
is on or near the boundary of the parameter space (Marschner, 2015).

Alternative methods for finding the MLE are provided in this package. For both GLMs and GAMs,
two approaches based on the EM algorithm can be used: a combinatorial EM (CEM) algorithm 
(Marschner, 2014) or an expanded EM algorithm. These accomodate the parameter constraints
and are more stable than iteratively reweighted least squares.

In a CEM algorithm, a collection of restricted parameter spaces is defined
which covers the full parameter space, and an EM algorithm is applied within each
restricted parameter space in order to find a collection of restricted maxima of
the log-likelihood function, from which can be obtained the global maximum over the
full parameter space. The methodology implemented for this algorithm is presented in 
Marschner and Gillett (2012) and Donoghoe and Marschner (2015).

In the expanded EM approach, additional parameters are added to the model,
and an EM algorithm finds the MLE of this overparameterised model by imposing
constraints on each individual parameter. This requires a single application of
the EM algorithm.

In each case, the EM algorithm may be accelerated by using the capabilities of the
\pkg{\link[turboEM:turboem]{turboEM}} package.

For GLMs, an adaptive barrier approach, which uses a constrained optimisation
algorithm, is also provided.
}
\author{
Mark W. Donoghoe \email{markdonoghoe@gmail.com}

Maintainer: Mark W. Donoghoe \email{markdonoghoe@gmail.com}
}
\references{
Donoghoe, M. W. and I. C. Marschner (2015). Flexible regression models for rate 
differences, risk differences and relative risks. \emph{International Journal of 
Biostatistics} 11(1): 91--108.

Donoghoe, M. W. and I. C. Marschner (2018). logbin: An R package for relative risk
regression using the log-binomial model. \emph{Journal of Statistical Software} 
86(9): 1--22.

Marschner, I. C. (2014). Combinatorial EM algorithms.
\emph{Statistics and Computing} 24(6): 921--940.

Marschner, I. C. (2015). Relative risk regression for binary outcomes: Methods
and recommendations. \emph{Australian & New Zealand Journal of Statistics}. In press.

Marschner, I. C. and A. C. Gillett (2012). Relative risk regression: Reliable
and flexible methods for log-binomial models. \emph{Biostatistics} 13(1): 179--192.
}
\keyword{ package }
\keyword{ regression }
\seealso{
\code{\link{glm}}, \pkg{\link[=glm2-package]{glm2}}, \pkg{\link[turboEM:turboem]{turboEM}}
}
\examples{
## For examples, see example(logbin) and example(logbin.smooth)
}
