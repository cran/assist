\name{sine4p}
\alias{sine4p}
\title{Calculate Reproducing Kernels for Periodic L-Splines with Period 1/2}
\description{
Return a matrix evaluating reproducing kernels for periodic L-splines at observed points.
}
\usage{
sine4p(s, t=s)
}
\arguments{
	\item{s}{a numeric vector.}
	\item{t}{an optional vector. Default is the same as s.}
}
\value{
a matrix with the numbers of row and column equal to the lengths of s and t respectively.
The [i, j] element is the reproducing kernel evaluated at (s[i], t[j]). }
\details{
The general formula of the reproducing kernel is provided in Gu (2001). The close form is not available, so an approximate based on the first 50 terms of the series
is used.}
\reference{
Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.

Gu, C. (2001). Smoothing Spline ANOVA Modes. Chapman and Hall.
}
\seealso{
\code{\link{cubic}}, \code{\link{lspline}}
}
\examples{
x<- seq(0, 1, len=100)
sine4p(x)
}
\keyword{file}


