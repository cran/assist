\name{Polynomial2}
\alias{linear2}
\alias{cubic2}
\alias{quintic2}
\alias{septic2}
\title{
Calculate Reproducing Kernels for Polynomial Splines on [0, T]
}
\description{
Return a matrix evaluating reproducing kernels for polynomial splines at observed points.
}
\usage{
linear2(s, t=s)
cubic2(s, t=s)
quintic2(s, t=s)
septic2(s, t=s)
}
\arguments{
   \item{s}{
   a vector of non-negative values, at which the kernels are evaluated.
   }
   \item{t}{
   an optional non-negative vector. Default is the same as s.
   }
}
\details{
The reproducing kernels implemented in these functions are based on Green functions. The domain is
[0, T], where T is a given positive number.
}
\value{
a matrix with the numbers of row and column equal to the length of s and t respectively.
The [i, j] element is the reproducing kernel of linear, cubic, quintic, or septic spline 
evaluated at (s[i], t[j]). 
}
\references{
Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.
}
\author{Chunlei Ke \email{chunlei\_ke@pstat.ucsb.edu} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{ssr}}, \code{\link{linear}}, \code{\link{cubic}}, 
\code{\link{quintic}}, \code{\link{septic}}
}
\examples{
x<- seq(0, 5, len=10)
linear2(x)
}
\keyword{file}

