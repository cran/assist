\name{sphere}
\alias{sphere}
\title{
Calculate Pseudo Reproducing Kernels for Spherical Splines
}
\description{
Return a matrix evaluating reproducing kernels for splines on a sphere.
}
\usage{
sphere(x, y=x, order=2)
}
\arguments{
\item{x}{
a matrix of two columns or a list of two components, representing observed 
latitude and longitude respectively.
}
\item{y}{
a matrix of two columns or a list of two components, representing 
latitude and longitude respectively. Default is the same as x.
}
\item{order}{
an optional integer sepcifying the order of the spherical spline. Available are
2, 3, 4, 5 and 6, with a default 2.
}}
\value{
a matrix with the numbers of row and column equal to the lengths of x and y respectively.
The [i, j] element is the reproducing kernel evaluated at \eqn{(x[i,], y[j,])} 
(or \eqn{((x[[1]][i], x[[2]][i]), (y[[1]][j], y[[2]][j]))} for lists). 
}
\details{
The kernel for sperical splines is a series inconvenient to compute. This pseudo kernel
is based on a topological equivalence as described in Wahba (1981), for which cases the
closed form can be derived.
}
\references{
Wahba, G. (1981). Spline Interprolation and Smoothing on the Sphere. SIAM J. Sci. Stat.Comput.,
Vol. 2, No. 1, March 1981.

Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.
}
\author{Chunlei Ke \email{chunlei_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{periodic}}
}
\examples{
\dontrun{
x<- seq(0, 2*pi, len=10)
y<- seq(-pi/2, pi/2, len=10)
s.ker<- sphere(cbind(x, y), order=3)
}
}
\keyword{file}
