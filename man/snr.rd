\name{snr}
\alias{snr}
\title{
Fit A Semi-parametric Nonlinear Regression Model
}
\description{
This generic function fits a Semi-parametric Nonlinear Regression Model as formulated in Ke (2000).
}
\usage{
snr(formula, func, params, data = sys.parent(), start, 
    spar = "v", verbose = FALSE, control = list(), correlation = NULL, 
    weights = NULL) 
}
\arguments{
\item{formula}{
a model formula, with the response on the left of a \eqn{\mbox{\textasciitilde}}{~} operator 
and on the right an expression representing the mean function 
with at least one unknown function appearing with a symbol, 
e.g. f. If "data" is present, all names except the nonparametric 
function(s) used in the formula should be defined as parameters 
or variables in the data frame.
}
\item{func}{
a list of spline formulae each specifying the spline components 
necessary to estimate each non-parametric function. On the left 
of a  \eqn{\mbox{\textasciitilde}}{~} operator of each component is the unknow function symbol(s) 
as well as its arguments, while the right side is a list of two 
components \code{nb}, an optional one-side formula for representing 
the null space's bases, and a required \code{rk} structure. \code{nb} and 
\code{rk} are similar to \code{formula} and \code{rk} in ssr. A missing \code{nb} 
denotes an empty null space.   
}
\item{params}{
a two-sided formula specifying models for the parameters. 
The syntax of \code{params} in \code{gnls} is adopted. See \code{gnls} for details.
}
\item{data}{
an optional data frame containing the variables named in model, 
params, correlation and weights. By default the variables are taken 
from the environment from which snr is called.	
}
\item{start}{
a numeric list with two components: "params=", a vector of the size of the length of the unknown parameters, 
providing inital values for the paramters, and "f=" a list of vectors or expressions which input inital values for the unknown functions. 
If the unknown functions appear linear in the model, the intial values then are not necessary. 
}
\item{spar}{
a character string specifying a method for choosing the smoothing parameter. "v", "m" and "u" represent  GCV, GML and
UBR respectively. Default is "v" for GCV.
}
\item{verbose}{
an optional logical numerical value. If \code{TRUE}, information on
the evolution of the iterative algorithm is printed. Default is
\code{TRUE}.
}
\item{control}{
an optional list of control parameters. See \code{snr.control} for details.
}
\item{correlation}{
an optional \code{corStruct} as in gnls. Default is NULL, corresponding to uncorrelation.
}
\item{weights }{
an optional \code{varFunc} structure as in \code{gnls}. Default is NULL, representing equal variances.
}}
\value{
An object of class \code{snr} is returned, representing a semi-parametric 
nonlinear regression fit. Generic functions such as print, summary, 
intervals and predict have methods to show the results of  the  fit.
}
\details{
A semi-parametric regression model is generalization of self-modeling 
regression, nonlinear regression and smoothing spline models, including 
as special cases (nonlinear) partial spline models, varying coefficients 
models, PP regression and some other popular models. 'snr' is 
implemented with an alternate iterative procedures with smoothing splines 
to estimate the unknown functions and general nonlinear regression to
estimate parameters.
}
\references{
Ke, C. (2000). Semi-parametric Nonlinear Regression and Mixed Effects 
Models. PhD thesis, University of California, Santa Barbara.

Pinheiro, J.C. and Bates, D. M. (2000). Mixed-Effects Models in S
and S-PLUS. Springer.

Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.
}
\author{Chunlei Ke \email{chunlei\_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}.}
\seealso{
\code{\link{intervals.snr}},  \code{\link{predict.snr}}, \code{\link{snr.control}}, 
\code{\link{gnls}}
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
}
\keyword{file}
