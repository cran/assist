\name{intervals.snr}
\alias{intervals.snr}
\title{
Calculate Predictions and Approximate Posterior Standard Deviations for Spline Estimates From a snr Object
}
\description{
Approximate posterior standard deviations are calculated for the spline estimate of
nonparametric functions from a snr object, based on which approximate Bayesian
confidence intervals may be constructed.
}
\usage{
intervals.snr(object, level=0.95, newdata=NULL, terms=<see below>, pstd=TRUE, level=0.95, ...)
}
\arguments{
\item{object}{
an object inheriting from class 'snr', representing a 
semi-parametric nonlinear regression model fit.
}
\item{level}{set as 0.95, unused currently}
\item{newdata}{
a data frame on which the fitted spline estimates are to be evaluated. 
Only those predictors, referred in 'func' of 'snr' fitting, have to
be present. The variable names of the data frame should correspond 
to the function(s)' arguments appearing in the opion func=  of snr.
Default is NULL, where predictions are made at the same values 
used to fit the object.
}
\item{terms }{
an optional named list of vectors or matrices containing 0's and 1's  collecting one or several combinations 
of the components of spline estimates in the fitted snr object. The length and names of the list shall match those of 
the unknown functions appearing in the 'snr' fit object. For the case of a single function, a vector of 0's 
and 1's can also be accepted. A value "1" at a particular position means that the component at 
that position is collected. Default is a vector of 1's, representing the overall fits of all unknown functions. 
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
a named list of objects of class "bCI" is returned, each component of which is a list of length 2. 
Within each component, the first element is a matrix which contains predictions for 
combinations specified by "terms", and the second element is a matrix which contains 
corresponding posterior standard deviations. 
}
\details{
The standard deviation returned is based on approximate Bayesian confidence 
intervals as formulated in Ke (2000). 
}
\references{
Ke, C. (2000). Semi-parametric Nonlinear Regression and Mixed Effects 
Models. PhD thesis, University of California, Santa Barbara.
}
\author{Chunlei Ke \email{chunlei\_ke@pstat.ucsb.edu} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{snr}}, \code{\link{plot.bCI}}, \code{\link{predict.ssr}}
}
\examples{
data(Arosa)
Arosa$csmonth <- (Arosa$month-0.5)/12
Arosa$csyear <- (Arosa$year-1)/45
# fit a simple sin-cos model
ozone.fit1 <- lm(thick~sin(2*pi*csmonth)+cos(2*pi*csmonth), data=Arosa)

# get start values
tmp <- atan(-ozone.fit1$coef[3]/ozone.fit1$coef[2])/(2*pi)
tmp <- log(tmp/(1-tmp))

# fit a SNR model
ozone.fit2 <- snr(thick~f1(csyear)+f2(csyear)*cos(2*pi*(csmonth+alogit(a))),
                   func=list(f1(x)+f2(x)~list(~x, cubic(x))), params=list(a~1),
                   data=Arosa, start=list(params=c(tmp)), spar="m")
p.ozone.fit2 <- intervals(ozone.fit2,data.frame(x=seq(0,1,len=100)),
            terms=list(f1=matrix(c(1,1,1,1,1,0,0,0,1),ncol=3,byrow=TRUE),
                       f2=matrix(c(1,1,1,1,1,0,0,0,1),ncol=3,byrow=TRUE)))

}
\keyword{file}

