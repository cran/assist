\name{Thin}
\alias{tp.pseudo}
\alias{tp}
\title{
Calculate Reproducing Kernels for Thin Plate Splines
}
\description{
Return a matrix evaluating reproducing kernels for thin plate splines at observed points.
}
\usage{
tp.pseudo(s, u=s, order=2)
tp(s, u=s, order=2)
}
\arguments{
\item{s}{
a list or matrix of observations. One component, if a list, and one column, if a matrix,
contains observations on one variable. If a list, all components must be of the same length. 
}
\item{u}{
a list or matrix of observations. If a list, all components must be of the same length. The number
of componets of the list, or the number of column of the matrix must be the same as that for s. 
Default is s.
}
\item{order}{
an optional integer specifying the order of the thin plate spline. Default is 2. Let d be the
dimension of s (and u). Then order must satisfy \eqn{2*order-d>0}.
}}
\value{
a matrix with the numbers of row and column equal to the common length of componets or 
the number of row of s and t respectively. The [i, j] element is the pseudo or true reproducing kernel
evaluated at the ith element of s and jth element of u.
}
\details{
The pseudo kernel, which is conditional definite positive instead of definite positive, is easy to
calculate, while the true reproducing kernel is complicated. Pseudo Kernels are enough to compute 
spline estimates, but to calcualte Bayesian confidnece intervals, the true kernel is required.
}
\references{
Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.

Gu, C. and Wahba, G (1993). Smoothing Spline ANOVA with component-wise Bayesian confidence intervals.
Journal of Computational and Graphical Statistics 55, 353--368.
}
\author{Chunlei Ke \email{chunlei\_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{ssr}}, \code{\link{cubic}}
}
\examples{
data(acid)
\dontrun{tp.pseudo(list(acid$x1, acid$x2))}
\dontrun{tp.pseud0(list(acid$x1, acid$x2), order=3)}
}
\keyword{file}
