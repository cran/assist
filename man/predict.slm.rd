\name{predict.slm}
\alias{predict.slm}
\title{
Predict Method for Semiparametric Linear Mixed Effects Model Fits 
}
\description{
Predicted Values on different levels of random effects with the spline fit
as part of fixed effects
}
\usage{
predict.slm(object, newdata=NULL, ...)
}
\arguments{
\item{object}{
an object inheriting from class \code{slm}, representing a semi-parametric linear
mixed effects model fit.
}
\item{newdata }{
a data frame containing the values at which predictions are required.
Only those predictors, referred to in the right side of the formula in
the object, need to be present by name in newdata. Default is NULL, where 
predictions are made at the same values used to  compute the object. 
}
\item{\dots}{other arguments, but currently unused.}
}
\value{
returned is a data.frame with  columns given by the predictions at different levels
and the grouping factors. Note that the smooth part of the spline fit is
regarded as fixed.
}
\references{
Wang, Y. (1998) Mixed Effects Smoothing Spline ANOVA. JRSS, Series B, 
60:159--174.

Pinherio, J. C. and Bates, D. M. (2000) Mixed-effects Models in S and S-Plus. Springer.
}
\author{Chunlei Ke \email{chunlei\_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}.}
\seealso{
\code{\link{slm}}
}

\examples{
data(dog)

dog.fit<-slm(y~group*time, rk=list(cubic(time), shrink1(group),
    rk.prod(kron(time-0.5),shrink1(group)),rk.prod(cubic(time), 
    shrink1(group))), random=list(dog=~1), data=dog)

predict(dog.fit)
}
\keyword{file}

