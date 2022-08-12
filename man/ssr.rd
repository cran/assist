\name{ssr}
\alias{ssr}
\title{Fit a General Smoothing Spline Regression Model}
\description{
Returns an object of class ssr which is a general/generalized/correlated smoothing spline fit. 
}
\usage{
ssr(formula, rk, data = list(), subset, weights = NULL, 
	correlation = NULL, family = "gaussian", scale = FALSE, 
	spar = "v", varht = NULL, limnla = c(-10, 3), control = list())
}
\arguments{
  \item{formula}{a \code{formula} object, with the response on the left of a \eqn{\mbox{\textasciitilde}}{~} operator, and 
the bases of the null space \eqn{H_0}, separated by + operators, on the right. 
Thus it specifies the parametric part of the model that contains functions 
which are not penalized. }
  \item{rk}{
a list of expressions specifying reproducing kernels \eqn{R^1},\dots,\eqn{R^p} for \eqn{H_1},\dots,\eqn{H_p}. 
For \eqn{p=1}, rk may be specified with given functions. Supported functions are: "linear", 
"cubic", "quintic", and "septic" for linear, cubic, quintic and septic polynomial 
splines with "linear2", "cubic2", "quintic2", and "septic2" for another construction;
"periodic" for periodic splines; "shrink0" and "shrink1" for Stein's shrink-toward-zero and 
shrink-toward-mean estimates; "tp" for thin-plate-splines; "lspline" for L-splines. 
For details on these kernels, see their help files. Users may also write their own functions.}
  \item{data}{a data frame containing the variables occurring in the formula and the \code{rk}. If this option is not specified, 
the variables should be on the search list. Missing values are not allowed. }
  \item{subset}{ an optional expression indicating which subset of the rows of the  data should be used in the fit.  
This can be a logical vector (which is replicated to have length equal to the number of observations), 
a numeric vector indicating which observation numbers are to be  included, or a character vector  
of the row names to be included.  All observations are included by default.}
  \item{weights}{ a vector or a matrix specifying known weights for weighted smoothing, or a varFunc structure 
specifying a variance function structure. Its length, if a vector, or its number of columns and rows, 
if a matrix, must be equal to the length of responses. See documentations of nlme for availabe 
varFunc structures. The default is that all weights are equal. }
  \item{correlation}{  a corStruct object describing the correlation structure for random errors. See documentations 
of corClasses for availabe correlation structures. Default is NULL for no correlation.}
  \item{family}{an optional string specifying the distribution family. Families supported are "binary", "binomial", 
"poisson", "gamma" and "gaussian" for Bernoulli, binomial, poisson, gamma and Gaussian distributions 
respectively. Default is "gaussian". }
  \item{scale}{an optional logical value. If `TRUE', all covariates appearing in "rk" will be scaled into 
interval [0, 1]. This transformation will affect predict.ssr. Default is FALSE. }
  \item{spar}{ a character string specifying a method for choosing the smoothing parameter. "v", "m" and "u" represent 
GCV, GML and UBR respectively. "u\eqn{\sim}{~}", only used for non-Gaussian families, specifies UBR with an estimated variance. 
Default is "v". }
  \item{varht}{ needed only when 'u' is chosen for 'method', which gives the fixed variance in calculation of the UBR function. 
Default is NULL for `family="gaussian"' and 1 for all other families. }
  \item{limnla}{a vector of length one or two, specifying a search range for log10(n*lambda), where lambda is the smoothing 
parameter and n is the sample size. If it is a single value, the smoothing parameter will be fixed at this value. 
This option is only applicable to spline smoothing with a single smoothing parameter. }
  \item{control}{a list of iteration and algorithmic constants. See ssr.control for details and default values. }
}
\details{
We adopt notations in Wahba (1990) for the general spline and smoothing spline ANOVA models. 
Specifically, the functional relationship between the predictor and independent variable is unknown 
and is assumed to be in a reproducing kernel Hilbert space H. H is decomposed into \eqn{H_0} and 
\eqn{H_1+...+H_p}, where the null space \eqn{H_0} is a finite dimensional space spanned by 
bases specified at the right side of \eqn{\mbox{\textasciitilde}}{~} in formula, and 
\eqn{H_1},\dots ,\eqn{H_p} are reproducing kernel Hilbert spaces with reproducing kernel specified in the list rk. 

The function is estimated from weighted penalized least square. ssr can be used to fit the general spline and smoothing spline ANOVA models (Wahba, 1990), generalized spline models (Wang, 1997) and correlated spline models (Wang, 1998). ssr can also fit partial spline model with additional parametric terms specified in the formula (Wahba, 1990).

ssr could be slow and memory intensive, especially for large sample size and/or when p is large. 
For fitting a cubic spline with CV or GCV estimate of the smoothing parameter, 
the S-Plus function \code{smooth.spline} is more efficient.

Components can be extracted using extractor functions predict, deviance, residuals, and summary. The output can be modified using update. 
}
\value{
an object of class \code{ssr} is  returned. See ssr.object for details. 

Note: output from earlier versions of \code{ssr} shows incorrect smoothing spline parameters for SSANOVA, which is corrected in this version.
}
\references{ 
Gu, C. (1989). RKPACK and its applications: Fitting smoothing spline models. Proceedings of the Statistical Computing Section, ASA, 42-51.

Gu, C. (2002). Smoothing Spline ANOVA. Spinger, New York.

Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.

Wang, Y. (1995). GRKPACK: Fitting Smoothing Spline ANOVA Models for Exponential Families. Communications in Statistics: Simulation and Computation, 24: 1037-1059.

Wang, Y. (1998) Smoothing Spline Models with Correlated Random Errors. JASA, 93:341-348.

Ke, C. and Wang, Y. (2002) ASSIST: A Suite of S-plus functions Implementing Spline smoothing Techniques. 
Available at: \url{https://yuedong.faculty.pstat.ucsb.edu/}}
\author{Yuedong Wang \email{yuedong@pstat.ucsb.edu} and Chunlei Ke \email{chunlei_ke@yahoo.com} }

\seealso{ \code{\link{deviance.ssr}}, \code{\link{hat.ssr}},  \code{\link{plot.ssr}}, \code{\link{ssr.control}},  
	\code{\link{predict.ssr}}, \code{\link{print.ssr}}, 
	\code{\link{ssr.object}}, \code{\link{summary.ssr}}, \code{\link{smooth.spline}}.}

\examples{
\dontrun{
library(MASS)
# fitting a cubic spline
fit1<- ssr(accel~times, data=mcycle, scale=T, rk=cubic(times))
summary(fit1)

# using GML to choose the smoothing parameter
fit2<- update(fit1, spar="m")

data(acid)
## fit an additive thin plate spline
acid.fit<- ssr( ph ~ t1 + x1 + x2, rk = list(tp(t1), tp(list(x1, x2))), 
        data = acid, spar = "m", scale = FALSE)
acid.fit
}
}
\keyword{file }



