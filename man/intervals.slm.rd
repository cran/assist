\name{intervals.slm}
\alias{intervals.slm}
\title{
Calculate Predictions and Posterior Standard Deviations of Spline Estimates From a slm Object
}
\description{
Provide a way  to calculate approximate posterior standard deviations and fitted 
values at any specified values for any combinations of elements of the spline 
estimate of nonparametric functions from a \code{slm} object, based on which 
approximate Bayesian confidence intervals may be constructed.
}
\usage{
intervals.slm(object, level=0.95, newdata=NULL, terms=<see below>, pstd=TRUE, level=0.95, ...)
}
\arguments{
\item{object}{
an object inheriting from class "slm", representing a semi-parametric nonlinear regression model fit.
}
\item{level}{set as 0.95, unused currently}
\item{newdata}{
an optional data frame on which the fitted spline estimate is to be evaluated. 
}
\item{terms }{
an optional vector of 0's and 1's collecting a combination of components, or a matrix of 0's and 1's 
collecting several combinations of components, in a fitted ssr object. All components include bases on 
the right side of \eqn{\mbox{\textasciitilde}}{~} in the formula and reproducing kernels in the rk list. Note that the first component 
is usually a constant function if it is not specifically excluded in the formula. A value "1" at a particular 
position means that the component at that position is collected. Default is a vector of 1's, 
representing the overall fit. 
}
\item{pstd}{
an optional logic value.
If TRUE (the default), the posterior standard deviations are calculated. 
Orelse, only the predictions are calculated.
Computation required for posterior standard deviations could be intensive. 
}
\item{level}{a numeric value set as 0.95.}
\item{\dots}{other arguments, currently unused.}
}
\value{
an object of class \code{bCI} is returned, which is a list of length 2. Its first element is a matrix which contains predictions for 
combinations specified by \code{terms}, and second element is a matrix which contains 
corresponding posterior standard deviations. 
}
\details{
The standard deviation returned is based on approximate Bayesian confidence intervals as formulated
in Wang (1998). 
}
\references{
Wang, Y. (1998). Mixed-effects smoothing spline ANOVA. Journal of the Royal Statistical Society, Series B 60, 159-174.
}
\author{Chunlei Ke \email{chunlei\_ke@pstat.ucsb.edu} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{slm}}, \code{\link{plot.bCI}}, \code{\link{predict.ssr}}
}
\examples{
data(dog)
# fit a SLM model with random effects for dogs
dog.fit<-slm(y~group*time, rk=list(cubic(time), shrink1(group),
    rk.prod(kron(time-0.5),shrink1(group)),rk.prod(cubic(time), 
    shrink1(group))), random=list(dog=~1), data=dog)

intervals(dog.fit)
}
\keyword{file}
