\name{predict.ssr}
\alias{predict.ssr}
\title{ 
Calculate Predictions and Posterior Standard Deviations for a ssr Object 
}
\description{ 
Provide a way to calculate predictions at any specified values for any combinations of elements in the fitted model. Posterior standard deviations may be used to construct Bayesian confidence intervals. 
}
\usage{ 
\method{predict}{ssr}(object, newdata=NULL, terms, pstd=TRUE, ...) 
}
\arguments{ 
 \item{object}{
  a fitted \code{ssr} object. 
  }
 \item{newdata}{
  an optional data frame containing the values at which predictions are required. Default is NULL, where predictions are made at the same values used to  compute the object. Note that if scale=T, the newdata is on the original scale before transformation. 
  }
  \item{terms }{
an optional vector of 0's and 1's collecting a combination of components, or a matrix of 0's and 1's collecting several combinations of components, in a fitted ssr object. All components include bases on the right side of \eqn{\mbox{\textasciitilde}}{~} in the formula and reproducing kernels in the rk list. Note that the first component is usually a constant function if it is not specifically excluded in the formula. A value "1" at a particular position means that the component at that position is collected. Default is a vector of 1's, representing the overall fit. 
  }
 \item{pstd}{
  an optional logic value. If TRUE (the default), the posterior standard deviations are calculated. Otherwise, only the predictions are calculated. Computation required for posterior standard deviations could be intensive. 
  }
\item{\dots}{other arguments, but currently unused.}
}
\details{
This function is a method for the generic function predict for class \code{ssr}. 
It can be used to construct Bayesian confidence intervals for any combinations 
of components in the fitted model. 
}
\value{ 
an object of class \code{bCI} is returned, which is a list of length 2. Its first element is a matrix which contains predictions for combinations specified by \code{terms}, and second element is a matrix which contains corresponding posterior standard deviations. 
}
\references{ 
Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.
}
\author{Chunlei Ke \email{chunlei\_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}.}
\seealso{
\code{\link{ssr}}, \code{\link{plot.bCI}}
}
\examples{
\dontrun{
data(acid)

# tp.pseudo calculates the pseudo kernel
acid.fit<- ssr( ph ~ t1 + x1 + x2, rk = list(tp.pseudo(t1), 
       tp.pseudo(list(x1, x2))), spar = "m", data=acid)

# extract the main effect of t1 
grid <- seq(min(acid$t1),max(acid$t1),length=100)
p <- predict(acid.fit,data.frame(t1=grid,x1=0,x2=0),
     terms=c(0,1,0,0,1,0))

# extract the main effect of (x1,x2) 
grid <- expand.grid(x1=seq(min(acid$x1),max(acid$x1),length=20),
     x2=seq(min(acid$x2),max(acid$x2),length=20))
p <- predict(acid.fit,data.frame(t1=0,x1=grid$x1,x2=grid$x2),
     terms=c(0,0,1,1,0,1),pstd=FALSE)
}
}
\keyword{file}
