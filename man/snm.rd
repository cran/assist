\name{snm}
\alias{snm}
\title{
Fit a Semi-parametric Nonlinear Mixed-effects Model
}
\description{
This generic function fits a semi-paramteric nonlinear mixed-effects model 
in the formulation described in Ke and Wang (2001). Current version only allows linear dependence on non-parametric 
functions.
}
\usage{
snm(formula, func, data=list(), fixed, random=fixed, 
groups, start, spar="v", verbose=FALSE, method="REML", control=NULL, 
correlation=NULL, weights=NULL)
}
\arguments{
\item{formula}{
a formula object, with the response on the left of a ~ operator, and an expression 
of variables, parameters and non-parametric functions on the right.
}
\item{func}{
a list of spline formulae each specifying the spline components necessary to 
estimate each non-parametric function. On the left of a \eqn{\mbox{\textasciitilde}}{~} operator of each component 
is the unknow function symbol(s) as well as its arguments, while the right side is a 
list of two components \code{nb}, an optional one-side formula for representing the null 
space's bases, and a required \code{rk} structure. \code{nb} and \code{rk} are similar to \code{formula} 
and \code{rk} in ssr. A missing \code{nb} denotes an empty null space. 
}
\item{fixed}{
a two-sided formula specifying models for the fixed effects.
The syntax of \code{fixed} in \code{nlme} is adopted.
}
\item{start}{
a numeric vector, the same length as the number of fixed effects, supplying starting
values for the fixed effects.
}
\item{spar}{
a character string specifying a method for choosing the smoothing parameter. "v", "m" and "u" represent  GCV, GML and
UBR respectively. Default is "v" for GCV.
}
\item{data}{
An optional data frame containing the variables appearing in \code{formula}
, \code{random}, \code{rk}, \code{correlation}, \code{weights}. By default, the variables 
are taken from the environment from which \code{snm} is called.
}
\item{random}{
an optional random effects structure specifying models for the random effects.
The same syntax of \code{random} in \code{nlme} is assumed.
}
\item{groups}{
an optional one-sided formula of the form ~g1 (single level) or ~g1/\dots/gQ 
(multiple levels of nesting), specifying the partitions of the data over 
which the random effects vary. g1,\dots,gQ must evaluate to factors in data. 
See nlme for details.	 
}
\item{verbose}{
an optional logical numerical value. If \code{TRUE}, information on
the evolution of the iterative algorithm is printed. Default is
\code{FALSE}.
}
\item{method}{
 a character string. If 'REML' the model is fit by maximizing the restricted 
log-likelihood. If 'ML' the log-likelihood is maximized. Default is 'REML. 
}
\item{control}{
a list of parameters to control the performance of the algorithm.
}
\item{correlation}{
an optional \code{corStruct} object describing the within-group correlation 
structure. See the documentation of corClasses for a description of the available corStruct classes. 
Default is NULL, corresponding to no within-in group correlations.
}
\item{weights}{
an optional \code{varFunc} object or one-sided formula describing the 
within-group heteroscedasticity structure. If given as a formula, 
it is used as the argument to \code{varFixed}, corresponding to fixed variance weights. 
See the documentation on varClasses for a description of the available varFunc 
classes. Defaults to NULL, corresponding to homoscesdatic within-group errors.
}}
\value{
an object of class \code{snm} is returned, representing a semi-parametric nonlinear
mixed effects model fit. Generic functions such as print, summary, predict and
intervals have methods to show the results of the fit.
}
\references{
Ke, C. and Wang, Y. (2001). Semi-parametric Nonlinear Mixed Effects Models and
Their Applications. JASA 96:1272-1298.

Pinheiro, J.C. and Bates, D. M. (2000). Mixed-Effects Models in S
and S-PLUS. Springer.
}
\author{Chunlei Ke \email{chunlei_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}.}
\seealso{
\code{\link{predict.snm}}, \code{\link{intervals.snm}}, \code{\link{snm.control}},
\code{\link{print.snm}},\code{\link{summary.snm}}
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
co2.fit2 <- snm(uptake~exp(a1)*f(exp(a2)*(conc-a3)),
                func=f(u)~list(~I(1-exp(-u))-1,lspline(u, type="exp")),
                fixed=list(a1~M-1,a3~1,a2~Type*Treatment),
                random=list(a1~1), group=~Plant, verbose=TRUE,
                start=co2.fit1$coe$fixed[c(2:4,9,5:8)], data=CO2)
}
}
\keyword{file}
