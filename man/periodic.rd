\name{periodic}
\alias{periodic}
\title{Calculate Reproducing Kernels for Periodic Polynomial Splines with Period 1}
\description{
Return a matrix evaluating reproducing kernels for periodic polynomial splines at observed points.
}
\usage{
periodic(s, t=s, order=2)
}
\arguments{
	\item{s}{a numeric vector.}
	\item{t}{an optional vector. Default is the same as s.}
	\item{order}{an optional integer sepcifying the order of the polynomial spline. Default is 2 for the 
periodic cubic spline.}
}
\value{
a matrix with the numbers of row and column equal to the lengths of s and t respectively.
The [i, j] element is the reproducing kernel evaluated at (s[i], t[j]). }
\details{
The general formula of the reproducing kernel is sum of an infinite series, which is approximated
by taking the first 50 terms. For the case of order=2, the close form is available and used.}
\references{
Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.

Gu, C. (2001). Smoothing Spline ANOVA Modes. Chapman and Hall.
}
\seealso{
\code{\link{cubic}}, \code{\link{lspline}}
}
\examples{
\dontrun{
x<- seq(0, 1, len=100)
periodic(x, order=3)
}
}
\keyword{file}


