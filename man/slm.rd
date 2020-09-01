\name{slm}
\alias{slm}
\title{
Fit a Semi-parametric Linear Mixed Effects Model
}
\description{
Returns an object of class \code{slm} that represents a
semi-parametric linear mixed effects model fit. 
}
\usage{
slm(formula, rk, data=list(), random, weights=NULL, 
correlation=NULL, control=list(apVar=FALSE))
}
\arguments{
\item{formula}{
a formula object, with the response on the left of a \eqn{\sim}{~} operator, and the bases 
of the null space \eqn{H_0} of the non-parametric function and other terms, separated by + operators, on the right.
}
\item{rk}{
a list of expressions that specify the reproducing kernels of the spline function(s), \eqn{R^1,\dots,R^p} for spaces \eqn{H_1,\dots,H_p}. See the help file of ssr for more details.
}
\item{data}{
An optional data frame containing the variables appearing in \code{formula}, \code{random}, \code{rk}, \code{correlation}, \code{weights}. 
By default, the variables are taken from the environment from which \code{slm} is called.
}
\item{random}{
A named list of formulae, lists of formulae, or pdMat objects, which defines
nested random effects structures. See help file of lme for more details.
}
\item{weights}{
An optional \code{varFun} object or one-sided formula describing the within-group heteroscedasticity stucture. 
See the help file of \code{lme} for more details.
}
\item{correlation}{
An optional \code{corStruct} object specifying the within-group correlation structure. See \code{lme} for more details.
}
\item{control}{
an optional list of any applicable control parameters from \code{lme}.
}}
\value{
An object of class \code{slm} is returned. Generic functions such as print, summary, predict and intervals have
methods to show the results of the fit.

Note: output from earlier versions of \code{slm} shows incorrect smoothing spline parameters for SSANOVA, which is corrected in this version.
}
\details{
This generic function fits a semi-parametric linear mixed effects model (or non-parametric mixed effects models) 
as described in Wang (1998), but allowing for general random and correlation structures. Because the connection
to a linear mixed effects model is adopted, only GML is available to choose smoothing parameters.
}
\references{
Wang, Y. (1998) Mixed Effects Smoothing Spline ANOVA. JRSS, Series B, 60:159--174.

Pinherio, J. C. and Bates, D. M. (2000) Mixed-effects Models in S and S-Plus. Springer.       
}
\author{Chunlei Ke \email{chunlei\_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}.}
\seealso{
\code{\link{ssr}}, \code{\link{predict.slm}}, \code{\link{intervals.slm}},
\code{\link{print.slm}},\code{\link{summary.slm}}
}
\examples{
\dontrun{
## SS ANOVA is used to model "time" and "group" 
## with random intercept for "dog".
data(dog)

dog.fit<- slm(y~group*time, rk=list(cubic(time), shrink1(group),
    	rk.prod(kron(time-0.5),shrink1(group)),rk.prod(cubic(time), 
    	shrink1(group))), random=list(dog=~1), data=dog)
}
}
\keyword{file}


