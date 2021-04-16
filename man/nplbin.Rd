\name{nplbin}
\alias{nplbin}
\title{
Non-Positive Log-Binomial Regression
}
\description{
Finds the maximum likelihood estimate of a log-link binomial GLM
using an EM algorithm, where each of the coefficients in the linear 
predictor is restricted to be non-positive.
}
\usage{
nplbin(y, x, offset, start, Amat = diag(ncol(x)), control = logbin.control(),
       accelerate = c("em", "squarem", "pem", "qn"),
       control.accelerate = list(list()))
}
\arguments{
  \item{y}{
	binomial response. May be a single column of 0/1
	or two columns, giving the number of successes and
	failures.
}
  \item{x}{
	non-negative covariate matrix.
}
  \item{offset}{
	non-positive additive offset vector. The default is a vector of zeros.
}
  \item{start}{
	starting values for the parameter estimates. All elements must be less than
	or equal to \code{-control$bound.tol}.
}
  \item{Amat}{
  matrix that parameter estimates are left-multiplied by before testing for
  convergence (e.g. to check reduced version of expanded parameter vector).
}
  \item{control}{
	a \code{\link{logbin.control}} object, which controls the fitting process.
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
}
\details{
This is a workhorse function for \code{\link{logbin}}, and runs the EM algorithm to find the 
constrained non-positive MLE associated with a log-link binomial GLM. See Marschner
and Gillett (2012) for full details.
}
\value{
A list containing the following components
\item{coefficients}{the constrained non-positive maximum likelihood estimate of the parameters.}
\item{residuals}{the residuals at the MLE, that is \code{y - fitted.values}}
\item{fitted.values}{the fitted mean values.}
\item{rank}{the number of parameters in the model (named "\code{rank}" for compatibility ---
we assume that models have full rank)}
\item{family}{included for compatibility --- will always be \code{\link{binomial}(log)}.}
\item{linear.predictors}{the linear fit on link scale.}
\item{deviance}{up to a constant, minus twice the maximised log-likelihood.}
\item{aic}{a version of Akaike's \emph{An Information Criterion},
minus twice the maximised log-likelihood plus twice the 
number of parameters.}
\item{aic.c}{a small-sample corrected
version of Akaike's \emph{An Information Criterion}
(Hurvich, Simonoff and Tsai, 1998).}
\item{null.deviance}{the deviance for the null model, comparable with \code{deviance}.
The null model will include the offset and an intercept.}
\item{iter}{the number of iterations of the EM algorithm used.}
\item{weights}{included for compatibility --- a vector of ones.}
\item{prior.weights}{the number of trials associated with each binomial response.}
\item{df.residual}{the residual degrees of freedom.}
\item{df.null}{the residual degrees of freedom for the null model.}
\item{y}{the \code{y} vector used.}
\item{converged}{logical. Did the EM algorithm converge?}
\item{boundary}{logical. Is the MLE on the boundary of the parameter
space --- i.e. are any of the \code{coefficients < control$bound.tol}?}
\item{loglik}{the maximised log-likelihood.}
\item{nn.design}{the non-negative \code{x} matrix used.}
}
\references{
Hurvich, C. M., J. S. Simonoff and C.-L. Tsai (1998). Smoothing parameter 
selection in non-parametric regression using an improved Akaike 
information criterion. \emph{Journal of the Royal Statistical Society: Series B 
(Statistical Methodology)} 60(2): 271--293.

Marschner, I. C. and A. C. Gillett (2012). Relative risk regression: reliable
and flexible methods for log-binomial models. \emph{Biostatistics} 13(1): 179--192.
}
\author{
Mark W. Donoghoe \email{markdonoghoe@gmail.com}.

This function is based on code from Marschner and Gillett (2012) written by Alexandra Gillett.
}
\keyword{ regression }
