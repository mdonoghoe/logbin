\name{predict.logbin.smooth}
\alias{predict.logbin.smooth}
\title{
Predict Method for logbin.smooth Fits
}
\description{
Obtains predictions from a fitted \code{\link{logbin.smooth}} object.
}
\usage{
\method{predict}{logbin.smooth}(object, newdata = NULL, type = c("link", "response", "terms"), 
        terms = NULL, na.action = na.pass, ...)
}
\arguments{
  \item{object}{
	a fitted object of class inheriting from \code{"logbin.smooth"}.
}
  \item{newdata}{
	optionally, a data frame in which to look for variables with which to
	predict. If omitted, the fitted linear predictors are used.
}
  \item{type}{
	the type of prediction required. The default is on the scale of the linear predictors; 
	the alternative \code{"response"} is on the scale of the response variable. 
	The \code{"terms"} option returns a matrix giving the fitted values of each term in the 
	model formula on the linear predictor scale.

	The value of this argument can be abbreviated.
}
  \item{terms}{
	with \code{type = "terms"} by default all terms are returned. 
	A character vector specifies which terms are to be returned.
}
  \item{na.action}{
	function determining what should be done with missing values in \code{newdata}. 
	The default is to predict \code{NA}.
}
  \item{\dots}{
	further arguments passed to or from other methods.
}
}
\details{
\code{predict.logbin.smooth} constructs the underlying basis functions for smooth variables
in \code{newdata} and runs \code{\link{predict.logbin}} to obtain predictions. Note that
if values of smooth covariates in \code{newdata} are outside the covariate space of 
\code{object}, an error will be returned.

If \code{newdata} is omitted, the predictions are based on the data used for the fit. 
In that case how cases with missing values in the original fit are treated is determined by the 
\code{na.action} argument of that fit. If \code{na.action = na.omit}, omitted cases 
will not appear in the residuals, whereas if \code{na.action = na.exclude} they will 
appear, with residual value \code{NA}. See also \code{\link{napredict}}.
}
\value{
A vector or matrix of predictions. For \code{type = "terms"}, this is a matrix with
a column per term, and may have an attribute \code{"constant"}.
}
\author{
Mark W. Donoghoe \email{markdonoghoe@gmail.com}
}
\note{
Variables are first looked for in \code{newdata} and then searched for in the usual 
way (which will include the environment of the formula used in the fit). A warning 
will be given if the variables found are not of the same length as those in 
\code{newdata} if it was supplied.
}
\seealso{
\code{\link{logbin.smooth}}, \code{\link{predict.logbin}}

\code{\link{predict.glm}} for the equivalent method for models fit using \code{\link{glm}}.
}
\examples{
## For an example, see example(logbin.smooth)
}
\keyword{ models }
\keyword{ regression }
\keyword{ smooth }
