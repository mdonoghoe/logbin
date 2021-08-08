\name{logbin}
\alias{logbin}
\title{
Log-Binomial Regression
}
\concept{Binomial regression}
\concept{CEM algorithm}
\description{
\code{logbin} fits relative risk (log-link) binomial
regression models.
}
\usage{
logbin(formula, mono = NULL, data, subset, na.action, start = NULL,
       offset, control = list(...), model = TRUE, 
       method = c("cem", "em", "glm", "glm2", "ab"),
       accelerate = c("em", "squarem", "pem", "qn"), 
       control.method = list(), warn = TRUE, ...)
}
\arguments{
  \item{formula}{
	an object of class \code{"\link{formula}"}
  (or one that can be coerced into that class): a symbolic
  description of the model to be fitted. The details of
  model specification are given under "Details".  Note that
  the model must contain an intercept, and 2nd-order terms
  (such as interactions) or above are currently not supported
  by the \code{"cem"} and \code{"em"} methods --- see "Note".
}
  \item{mono}{
	a vector indicating which terms in
  \code{formula} should be restricted to have a
  monotonically non-decreasing relationship with the
  outcome. May be specified as names or indices of the
  terms.
  
  \code{method = "glm"} and \code{"glm2"} cannot impose
  monotonicity constraints, and they are not currently
  supported for \code{method = "ab"}.
}
  \item{data}{
	an optional data frame, list or environment
  (or object coercible by \code{\link{as.data.frame}} to a
  data frame) containing the variables in the model. If not
  found in \code{data}, the variables are taken from
  \code{\link{environment}(formula)}, typically the environment
  from which \code{logbin} is called.
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
  \item{start}{
	starting values for the parameters in the
  linear predictor.
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
  
  With \code{method = "cem"}, \code{epsilon}
  should be smaller than \code{bound.tol}.
}
  \item{model}{
	a logical value indicating whether the \emph{model frame} 
	should be included as a component of the returned value.
}
  \item{method}{
    a character string that determines which algorithm to use
  to find the MLE. The main purpose of \code{logbin} is the
  implementation of stable EM-type algorithms: \code{"cem"} for 
  the combinatorial EM algorithm, which cycles through a sequence 
  of constrained parameter spaces, or \code{"em"} for a single EM 
  algorithm based on an overparameterised model.
  
  \code{"ab"} implements an adaptive barrier method, using the
  \code{\link{constrOptim}} function.
  
  \code{"\link{glm}"} or \code{"\link[glm2]{glm2}"} may be used
  to compare the results from the usual IWLS algorithms on the same model.
}
  \item{accelerate}{
    for the \code{"cem"} and \code{"em"} methods, a character string that determines the acceleration
  algorithm to be used, (partially) matching one of \code{"em"} (no acceleration --- the default), 
  \code{"squarem"}, \code{"pem"} or \code{"qn"}. See \code{\link[turboEM]{turboem}}
  for further details. Note that \code{"decme"} is not permitted.
}
  \item{control.method}{
    a list of control parameters for the fitting algorithm.
    
  This is passed to the \code{control.method} argument of \code{\link[turboEM]{turboem}}
  if \code{method = "cem"} or \code{"em"}.
  
  If \code{method = "ab"}, this is passed to the \code{control} argument of
  \code{\link{constrOptim}} (and hence to \code{\link{optim}} --- see this documentation
  for full details). Note that the \code{trace} and \code{maxit} elements are 
  ignored and the equivalent items from the supplied \code{logbin.control} argument
  are used instead. May also contain element \code{method} (default \code{"BFGS"}), 
  which is passed to the \code{method} argument of \code{\link{constrOptim}}.
  
  If any items are not specified, the defaults are used.
}
  \item{warn}{
	a logical indicating whether or not warnings should be
  provided for non-convergence or boundary values.
}
  \item{\dots}{
	arguments to be used to form the default
  \code{control} argument if it is not supplied directly.
}
}
\details{
\code{logbin} fits a generalised linear model (GLM) with a
binomial error distribution and log link
function. Predictors are assumed to be continuous, unless
they are of class \code{\link{factor}}, or are character or
logical (in which case they are converted to
\code{factor}s). Specifying a predictor as monotonic using
the \code{mono} argument means that for continuous terms,
the associated coefficient will be restricted to be
non-negative, and for categorical terms, the coefficients
will be non-decreasing in the order of the factor
\code{\link{levels}}. This allows semi-parametric monotonic regression
functions, in the form of unsmoothed step-functions. For 
smooth regression functions see \code{\link{logbin.smooth}}.

As well as allowing monotonicity constraints, the function 
is useful when a standard GLM routine, such as
\code{\link{glm}}, fails to converge with a log-link
binomial model. For convenience in comparing convergence on
the same model, \code{logbin} can be used
as a wrapper function to \code{\link{glm}} and \code{\link[glm2]{glm2}}
through the \code{method} argument. 

If \code{\link{glm}} does achieve successful convergence,
and \code{logbin} converges to an interior point, then the two
results will be identical. However, as illustrated in one of 
the examples below, \code{\link{glm}} may still experience convergence 
problems even when \code{logbin} converges to an interior point. 
Note that if \code{logbin} converges to a boundary point, then it 
may differ slightly from \code{\link{glm}} even if \code{\link{glm}} successfully 
converges, because of differences in the definition of the parameter 
space. \code{logbin} produces valid fitted values for covariate 
values within the Cartesian product of the observed range of covariate 
values, whereas \code{\link{glm}} produces valid fitted values just 
for the observed covariate combinations (assuming it successfully 
converges). This issue is only relevant when \code{logbin} 
converges to a boundary point. The adaptive barrier approach defines 
the parameter space in the same way as \code{\link{glm}}, so the
same comments apply when comparing its results to those from
\code{method = "cem"} or \code{"em"}.

The main computational method is an EM-type algorithm which accommodates 
the parameter contraints in the model and is more stable than iteratively
reweighted least squares. This is done in one of two ways,
depending on the choice of the \code{method} argument.

\code{method = "cem"} implements a CEM algorithm (Marschner, 2014),
in which a collection of restricted parameter spaces is defined 
that covers the full parameter space, and an EM algorithm is applied within each
restricted parameter space in order to find a collection of
restricted maxima of the log-likelihood function, from
which can be obtained the global maximum over the full
parameter space. See Marschner and Gillett (2012) for further
details.

\code{method = "em"} implements a single EM algorithm
on an overparameterised model, and the MLE of this model
is transformed back to the original parameter space.

Acceleration of the EM algorithm in either case can be
achieved through the methods of the \code{\link[turboEM]{turboem}}
package, specified through the \code{accelerate} argument. However,
note that these methods do not have the guaranteed convergence of
the standard EM algorithm, particularly when the MLE is on the
boundary of its (possibly constrained) parameter space.

Alternatively, an adaptive barrier method can be used by specifying
\code{method = "ab"}, which maximises the likelihood subject to
constraints on the fitted values.
}
\value{
\code{logbin} returns an object of class \code{"logbin"},
which inherits from classes \code{"glm"} and \code{"lm"}.
The function \code{\link{summary.logbin}} can be used 
to obtain or print a summary of the results.

The generic accessor functions \code{\link{coefficients}},
\code{\link{fitted.values}} and \code{\link{residuals}} can be used to
extract various useful features of the value returned by
\code{logbin}. Note that \code{\link{effects}} will not work.

An object of class \code{"logbin"} is a list containing the
same components as an object of class \code{"glm"} (see the
"Value" section of \code{\link{glm}}). It also includes: 
\item{loglik}{the maximised log-likelihood.}
\item{aic.c}{a small-sample corrected
version of Akaike's \emph{An Information Criterion}
(Hurvich, Simonoff and Tsai, 1998). This is used by
\code{\link{logbin.smooth}} to choose the optimal number of
knots for smooth terms.} 
\item{xminmax}{the minimum and maximum observed values
for each of the continuous covariates, to help define
the covariate space of the model.}

As well as:
\item{np.coefficients}{estimated coefficients associated
with the non-positive parameterisation corresponding to the
MLE.} 
\item{nn.x}{non-negative model matrix associated with
\code{np.coefficients}.}
\item{coefhist}{(if \code{control$coeftrace = TRUE}), a matrix or
list of matrices containing the coefficient estimates after each
EM iteration.}
}
\references{
Hurvich, C. M., J. S. Simonoff and C.-L. Tsai (1998). Smoothing parameter 
selection in non-parametric regression using an improved Akaike 
information criterion. \emph{Journal of the Royal Statistical Society: Series B 
(Statistical Methodology)} 60(2): 271--293.

Donoghoe, M. W. and I. C. Marschner (2018). logbin: An R package for relative risk
regression using the log-binomial model. \emph{Journal of Statistical Software} 
86(9): 1--22.

Marschner, I. C. (2014). Combinatorial EM algorithms.
\emph{Statistics and Computing} 24(6): 921--940.

Marschner, I. C. and A. C. Gillett (2012). Relative risk regression: reliable
and flexible methods for log-binomial models. \emph{Biostatistics} 13(1): 179--192.
}
\author{
Mark W. Donoghoe \email{markdonoghoe@gmail.com}
}
\note{
Due to the way in which the covariate space is defined in the CEM algorithm,
models that include terms that are functionally dependent on one another
--- such as interactions and polynomials --- may give unexpected
results. Categorical covariates should always be entered directly
as factors rather than dummy variables. 2-way interactions between 
factors can be included by calculating a new factor term that 
has levels corresponding to all possible combinations of the factor
levels (see the Example). Non-linear relationships can be included
by using \code{\link{logbin.smooth}}.
}
\seealso{
\code{\link{logbin.smooth}} for semi-parametric models

\code{\link[turboEM]{turboem}} for acceleration methods

\code{\link{constrOptim}} for the adaptive barrier approach.
}
\examples{
require(glm2)
data(heart)

#======================================================
#  Model with periodic non-convergence when glm is used
#======================================================

start.p <- sum(heart$Deaths) / sum(heart$Patients)

fit.glm <- glm(cbind(Deaths, Patients-Deaths) ~ factor(AgeGroup) + factor(Severity) +
  factor(Delay) + factor(Region), family = binomial(log), 
  start = c(log(start.p), rep(c(0.2, 0.4), 4)), data = heart,
  trace = TRUE, maxit = 100)

fit.logbin <- logbin(formula(fit.glm), data = heart, 
                     start = c(log(start.p), rep(c(0.2, 0.4), 4)),
                     trace = 1)
summary(fit.logbin)

# Speed up convergence by using single EM algorithm
fit.logbin.em <- update(fit.logbin, method = "em")

# Speed up convergence by using acceleration methods
fit.logbin.acc <- update(fit.logbin, accelerate = "squarem")
fit.logbin.em.acc <- update(fit.logbin.em, accelerate = "squarem")

\donttest{
#=============================
#  Model with interaction term
#=============================

heart$AgeSev <- 10 * heart$AgeGroup + heart$Severity

fit.logbin.int <- logbin(cbind(Deaths, Patients-Deaths) ~ factor(AgeSev) +
  factor(Delay) + factor(Region), data = heart, trace = 1, maxit = 100000)
  
summary(fit.logbin.int)
vcov(fit.logbin.int)
confint(fit.logbin.int)
summary(predict(fit.logbin.int, type = "response"))
anova(fit.logbin, fit.logbin.int, test = "Chisq")
}
}
\keyword{ models }
\keyword{ regression }
