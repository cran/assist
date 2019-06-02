\name{intervals.nnr}
\alias{intervals.nnr}
\title{
Calculate Predictions and Approximate Posterior Standard Deviations for Spline Estimates From a nnr Object
}
\description{
Approximate posterior standard deviations are calculated for the spline estimate of
nonparametric functions from a \code{nnr} object, based on which approximate Bayesian
confidence intervals may be constructed.
}
\usage{
\method{intervals}{nnr}(object,level=0.95, newdata=NULL, terms, pstd=TRUE, ...)
}
\arguments{
\item{object}{
an object inheriting from class \code{nnr}, representing a 
nonlinear nonparametric regression model fit.
}
\item{newdata}{
a data frame on which the fitted spline estimates are to be evaluated. 
Only those predictors, referred in \code{func} of \code{nnr} fitting, have to
be present. The variable names of the data frame should correspond 
to the function(s)' arguments appearing in the opion func=  of nnr.
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
an object of class \code{bCI} is returned, which is a list of length 2. 
Its first element is a matrix which contains predictions for 
combinations specified by \code{terms}, and second element is a matrix which contains 
corresponding posterior standard deviations. 
}
\details{
The standard deviation returned is based on approximate Bayesian confidence 
intervals as formulated in Ke and Wang (2002). 
}
\references{
Ke, C. and Wang, Y. (2002). Nonlinear Nonparametric Regression Models. Submitted.
}
\author{Chunlei Ke \email{chunlei\_ke@pstat.ucsb.edu} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{nnr}}, \code{\link{plot.bCI}}
}
\examples{
\dontrun{
## fit a generalized varying coefficient models
data(Arosa)
Arosa$csmonth <- (Arosa$month-0.5)/12
Arosa$csyear <- (Arosa$year-1)/45
ozone.fit <- nnr(thick~f1(csyear)+exp(f2(csyear))*f3(csmonth),
        func=list(f1(x)~list(~I(x-.5),cubic(x)), f2(x)~list(~I(x-.5)-1,cubic(x)),
        f3(x)~list(~sin(2*pi*x)+cos(2*pi*x)-1,lspline(x,type="sine0"))),
        data=Arosa[Arosa$year\%\%2==1,], spar="m", start=list(f1=mean(thick),f2=0,f3=sin(csmonth)),
	control=list(backfit=1))

x <- seq(0,1,len=50)
u <- seq(0,1,len=50)

## calculate Bayesian confidence limits for all components of all functions
p.ozone.fit <- intervals(ozone.fit, newdata=list(csyear=x,csmonth=u),
                 terms=list(f1=matrix(c(1,1,1,1,1,0,0,0,1),nrow=3,byrow=TRUE),
	                    f2=matrix(c(1,1,1,0,0,1),nrow=3,byrow=TRUE),
                            f3=matrix(c(1,1,1,1,1,0,0,0,1),nrow=3,byrow=TRUE)))	
plot(p.ozone.fit, x.val=x)
}
}
\keyword{file}

