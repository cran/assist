\name{snr}
\alias{snr}
\title{
Fit A Semi-parametric Nonlinear Regression Model
}
\description{
This generic function fits a Semi-parametric Nonlinear Regression Model as formulated in Ke (2000).
}
\usage{
snr(formula, func, params, data, start, 
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
\code{\link{intervals.snr}},  \code{\link{predict.snr}}, \code{\link{snr.control}}
}
\examples{
\dontrun{
data(CO2)
options(contrasts=rep("contr.treatment", 2))    
co2.fit1 <- nlme(uptake~exp(a1)*(1-exp(-exp(a2)*(conc-a3))), 
                 fixed=list(a1+a2~Type*Treatment,a3~1), 
                 random=a1~1, groups=~Plant, 
                 start=c(log(30),0,0,0,log(0.01),0,0,0,50),
                 data=CO2)

M <- model.matrix(~Type*Treatment, data=CO2)[,-1]

## fit a SNR model
co2.fit2 <- snr(uptake~exp(a1)*f(exp(a2)*(conc-a3)),
                func=f(u)~list(~I(1-exp(-u))-1,lspline(u, type="exp")),
                params=list(a1~M-1, a3~1, a2~Type*Treatment),
                start=list(params=co2.fit1$coe$fixed[c(2:4,9,5:8)]), data=CO2)
}
}
\keyword{file}
