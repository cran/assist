\name{hat.ssr}
\alias{hat.ssr}
\title{
Extract the Hat Matrix from a ssr Object}
}
\description{
Calculate the hat matrix for a \code{ssr} object.
}
\usage{ 
hat.ssr(ssr.obj)
}
\arguments{
   \item{ssr.obj}{
   a fitted ssr object.
   }
}
\details{
The hat matrix may be used for diagnosis. Note that the full name hat.ssr shoud be used since the function hat already exist.
}
 
\value{
returned is the hat (influence, smoother) matrix.
}
\references{ 
Eubank, R. L. (1984). The Hat Matrix for Smoothing Splines. Statistics and Probability Letters, 2:9-14.

Eubank, R. L. (1985). Diagnostics for Smoothing Splines. Journal of the Royal Statistical Society B. 47: 332-341.

Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.
}
\author{Chunlei Ke \email{chunlei\_ke@pstat.ucsb.edu} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{ssr}}
}
\examples{
\dontrun{library(MASS)}
\dontrun{fit1<- ssr(accel~times, data=mcycle, scale=T, rk=cubic(times))}
\dontrun{h <- hat.ssr(fit1)}
}
\keyword{file}

