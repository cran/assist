\name{plot.bCI}
\alias{plot.bCI}
\title{
Bayesian Confidence Interval Plot of a Smoothing Spline Fit
}
\description{
Create trellis plots of a nonparametric function fit together
with its (approximate) 95\% Bayesian confidence intervals from 
a ssr/slm/snr/snm object.
}
\usage{
\method{plot}{bCI}(x, x.val=NULL, type.name=NULL, \dots)
}
\arguments{
\item{x}{
an object of class "bCI" containing point evaluation of the unknown
function and/or corresponding posterior standard devaitions.
}
\item{x.val}{
an optional vector representing values of argument based on which the 
function is to evaluate.
}
\item{type.name}{
an optional character vector specifying the names of fits.
}
\item{\dots}{
options suitable for xyplot.
}}
\details{
This function is to visualize a spline fit by use of trellis graphic facility 
with Bayesian confidence intervals superposed. Multi-panel plots, based on xyplot,
are suitable for SS ANOVA decomposition of a spline estimate. 
}
\author{Chunlei Ke \email{chunlei\_ke@pstat.ucsb.edu} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{predict.ssr}}, \code{\link{intervals.slm}}, 
\code{\link{intervals.snr}}, \code{\link{intervals.snm}}
}
\examples{
\dontrun{
x<- seq(0, 1, len=100)
y<- 2*sin(2*pi*x)+rnorm(x)*0.5

fit<- ssr(y~x, cubic(x))
p.fit<- predict(fit)
plot(p.fit)
plot(p.fit,type.name="fit")
}
}
\keyword{file}
