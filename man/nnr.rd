\name{nnr}
\alias{nnr}
\title{
Nonlinear Non-parametric Regression
}
\description{
Fit a nonlinear nonparametric regression models with spline smoothing based on extended Gauss-Newton/Newton-Raphson and backfitting.
}
\usage{
nnr(formula, func, spar="v", data=list(),
    start=list(),verbose=FALSE,  control=list())
}
\arguments{
\item{formula}{
a model formula, with the response on the left of a \eqn{\mbox{\textasciitilde}}{~} operator and on the right an expression representing 
the mean function with a nonparametric function appearing with a symbol, e.g. f. 
}
\item{func}{
a required formula specifying the spline components necessary to estimate the non-parametric function. 
On the left of a \eqn{\mbox{\textasciitilde}}{~} operator is the unknow function symbol as well as its arguments, while the right side 
is a list of two components, an optional \code{nb} and a required \code{rk}. \code{nb} and \code{rk} are 
similar to \code{formula} and \code{rk} in \code{ssr}. A missing \code{nb} denotes an empty null space.   
}
\item{spar}{
a character string specifying a method for choosing the smoothing parameter. "v", "m" and "u" represent  GCV, GML and
UBR respectively. Default is "v" for GCV.
}
\item{data}{
an optional data frame.
}
\item{start}{
a list of vectors or expressions which input inital values for the unknown functions. If expressions,
the argument(s) inside should be the same as in \code{func}. The length of \code{start} should be the same as 
the number of unknown functions. If named, the names of the list should match those in "func". If not named, the order 
of the list is taken as that appearing in "func".
}
\item{verbose}{
an optional logical numerical value. If \code{TRUE}, information on
the evolution of the iterative algorithm is printed. Default is \code{FALSE}.
}
\item{control}{
an optional list of control values to be used.  See nnr.control for details.
}}
\value{
an object of class \code{nnr} is returned, containing fitted values, fitted function values as well as 
other information used to assess the estimate.
}
\details{
A nonlinear nonparametric model is fitted using the algorithms developed in Ke and Wang (2002).
}
\references{
Ke, C. and Wang, Y. (2002). Nonlinear Nonparametric Regression Models. Submitted.
}
\author{Chunlei Ke \email{chunlei_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}.}
\seealso{
\code{\link{nnr.control}}, \code{\link{ssr}}, \code{\link{print.nnr}}, \code{\link{summary.nnr}}, \code{\link{intervals.nnr}}
}
\examples{
\dontrun{
x<- 1:100/100
y<- exp(sin(2*pi*x))+0.3*rnorm(x)
fit<- nnr(y~exp(f(x)), func=list(f(u)~list(~u, cubic(u))), start=list(0))

## fit a generalized varying coefficient models
data(Arosa)
Arosa$csmonth <- (Arosa$month-0.5)/12
Arosa$csyear <- (Arosa$year-1)/45
ozone.vc.fit <- nnr(thick~f1(csyear)+exp(f2(csyear))*f3(csmonth),
        func=list(f1(x)~list(~I(x-.5),cubic(x)), f2(x)~list(~I(x-.5)-1,cubic(x)),
        f3(x)~list(~sin(2*pi*x)+cos(2*pi*x)-1,lspline(x,type="sine0"))),
        data=Arosa[Arosa$year\%\%2==1,], spar="m", start=list(f1=mean(thick),f2=0,f3=sin(csmonth)),
        control=list(backfit=1))
}
}
\keyword{file}






