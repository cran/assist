\name{predict.snm}
\alias{predict.snm}
\title{
Predictions from a Semiparametric Nonlinear Mixed Effects Model Fit
}
\description{
The predictions are obtained on a semiparametric nonlinear mixed effects model object 
by replacing the unknown functuons and the unknown parameters with their estimates. 
Of note, only a population level of predictions is available. 
}
\usage{
\method{predict}{snm}(object, newdata, ...)
}
\arguments{
\item{object}{
a fitted \code{snm} object.
}
\item{newdata}{
a data frame containing the values at which predictions are required. 
Default are data used to fit the object.
}
\item{\dots}{other arguments, but currently unused.}
}
\value{
a vector of prediction values, obtained by evaluating the model in the frame \code{newdata}
}
\details{
This function is a method for the generic function predict for class \code{snm}. 
}
\references{
Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.

Ke, C. and Wang, Y. (2001). Semi-parametric Nonlinear Mixed Effects Models and
Their Applications. JASA.
}
\author{Chunlei Ke \email{chunlei\_ke@pstat.ucsb.edu} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{snm}}, \code{\link{predict}}
}
\keyword{file}


