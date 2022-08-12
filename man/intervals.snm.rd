\name{intervals.snm}
\alias{intervals.snm}
\title{
Calculate Predictions and Approximate Posterior Standard Deviations for Spline Estimate From a snm Object
}
\description{
Provide a way  to calculate approximate posterior standard deviations and fitted 
values at any specified values for any combinations of elements of the spline 
estimate of nonparametric functions from a snm object, based on which 
approximate Bayesian confidence intervals may be constructed.
}
\usage{
\method{intervals}{snm}(object,level=0.95,newdata=NULL, terms, pstd=TRUE, ...)
}
\arguments{
\item{object}{
an object inheriting from class snm, representing a semi-parametric 
nonlinear mixed effects model fit.
}
\item{newdata}{
a data frame on which the fitted spline estimates are to be evaluated. 
Only those predictors, referred in 'func' of 'snm' fitting, have to
be present. The variable names of the data frame should correspond 
to the function(s)' arguments appearing in the opion func=  of snm.
Default is NULL, where predictions are made at the same values 
used to fit the object.
}
\item{terms }{
an optional vector of 0's and 1's collecting a combination of 
components, or a matrix of 0's and 1's collecting several combinations 
of components of spline estimates in a fitted snm object. Note that 
in the cases of multiple functions, the order of all componets is 
collection of base functions for all functions followed by RK's. 
A value "1" at a particular position means that the component at 
that position is collected. Default is a vector of 1's, 
representing the overall fit. 
}
\item{pstd}{
an optional logic value. If TRUE (the default), approximate posterior 
standard deviations are calculated. Orelse, only the predictions 
are calculated. Computation required for posterior standard deviations 
could be intensive. 
}
\item{level}{a numeric value set as 0.95.}
\item{\dots}{other arguments, currently unused.}
}
\value{
an object of class \code{bCI} is returned, which is a list of length 2. 
Its first element is a matrix which contains predictions for 
combinations specified by "terms", and second element is a matrix 
which contains corresponding posterior standard deviations. 
}
\details{
The standard deviation returned is based on approximate Bayesian 
confidence intervals as formulated in Ke and Wang (2001). 
}
\references{
Ke, C. and Wang, Y. (2001). Semi-parametric Nonlinear Mixed Effects 
Models and Their Applications. JASA 96:1272-1298.
}
\author{Chunlei Ke \email{chunlei_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}.}
\seealso{
\code{\link{snm}}, \code{\link{plot.bCI}}, \code{\link{predict.ssr}}
}
\examples{
\dontrun{
data(horm.cort)

## extract normal dubjects
cort.nor<- horm.cort[horm.cort$type=="normal",]

## fit a self-modelling model with random effects
cort.fit<- snm(conc~b1+exp(b2)*f(time-alogit(b3)), 
  func=f(u)~list(periodic(u)), fixed=list(b1~1), 
  random=pdDiag(b1+b2+b3~1), data=cort.nor, 
  groups= ~ID,start=mean(cort.nor$conc))

## note the variable name of newdata
intervals(cort.fit, newdata=data.frame(u=seq(0,1,len=50)))
}
}
\keyword{file}

