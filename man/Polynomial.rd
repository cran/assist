\name{Polynomial}
\alias{linear}
\alias{cubic}
\alias{quintic}
\alias{septic}
\title{
Calculate Reproducing Kernels for Polynomial Splines on [0, 1]
}
\description{
Return a matrix evaluating reproducing kernels for polynomial splines at observed points.
}
\usage{
linear(s, t=s)
cubic(s, t=s)
quintic(s, t=s)
septic(s, t=s)
}
\arguments{
 \item{s}{
   a vector of values in [0, 1], at which the kernels are evaluated.
 }
 \item{t}{
   an optional vector in [0, 1]. Default is the same as s.
 }
}

\details{
The reproducing kernels implemented in these functions are based on Bernoulli functions 
with domain [0, 1].
}
\value{
a matrix with the numbers of row and column equal to the lengths of s and t respectively.
The [i, j] element is the reproducing kernel of linear, cubic, quintic, or septic spline 
evaluated at (s[i], t[j]). 
}
\references{
Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.
}
\author{Chunlei Ke \email{chunlei\_ke@pstat.ucsb.edu} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{ssr}}, \code{\link{linear2}}, \code{\link{cubic2}}, 
\code{\link{quintic2}}, \code{\link{septic2}}
}
\examples{
\dontrun{
x<-seq(0, 1, len=10)
cubic(x)
}
}
\keyword{file}
