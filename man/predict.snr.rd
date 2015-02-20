\name{predict.snr}
\alias{predict.snr}
\title{
Predict Method from a Semiparametric Nonlinear Regression Model Fit
}
\description{
The predictions on a semiparametric nonlinear regression model object are obtained by 
substituting the unknwon functions together with unknown parameters with their estimates 
and evaluating the regression functional based on provided or default covariate values.
}
\usage{
\method{predict}{snr}(object, newdata, ...)
}
\arguments{
\item{object}{
a fitted \code{snr} object.
}
\item{newdata}{
a data frame containing the values at which predictions are required. 
Default are NULL, where data used to produce the fit are to be taken.
}
\item{\dots}{other arguments, but currently unused.}
}
\value{
a vector of prediction values, obtained by evaluating the model in the frame \code{newdata}.
}
\details{
This function is a method for the generic function predict for class \code{snr} 
}
\references{
Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.


Ke, C. (2000). Semi-parametric Nonlinear Regression and Mixed Effects 
Models. PhD thesis, University of California, Santa Barbara.
}
\author{Chunlei Ke \email{chunlei\_ke@pstat.ucsb.edu} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{snr}}
}
\keyword{file}


