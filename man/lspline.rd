\name{lspline}
\alias{lspline}
\title{
Calculate Reproducing Kernels for Some L-splines
}
\description{
Return a matrix evaluating reproducing kernels for some L-splines at observed points.
}
\usage{
lspline(x,y=x, type="exp", \dots)
}
\arguments{
\item{x}{
a numeric vector on which reproducing kerenls are evaluated.
}
\item{y}{
an optional vector, specifying the second argument of reproducing kernels. Default is \code{x}.
}
\item{type}{
a string indicating the type of L-splines. Available options 
are "exp", "logit","sine", "sine1", and "linSinCos". Default is "exp".
}
\item{\dots}{
other arguments needed.}
}
\value{
a matrix with the numbers of row and column equal to the lengths of x and y respectively.
The [i, j] element is the reproducing kernel evaluated at (x[i], y[j]). 
}
\details{
Denote L as the differential oprator, \eqn{H_0} as the null (kernel) space. The available kernels
correspond to the following L:
\itemize{
    \item exp: \eqn{L=rD+D^2}, \eqn{H_0=span\{1,exp(-rx)\}}. \eqn{r>0}, default to be 1;\cr
    \item logit: \eqn{L=D-1/(1+e^t)}, \eqn{H_0=span\{e^t/(1+e^t)\}};\cr
    \item sine0: \eqn{L=D^2+(2\pi)^2}, \eqn{H_0=span\{sin(2\pi x),cos(2\pi x)\}};\cr
    \item sine1: \eqn{L=D(D^2+(2\pi)^2)}, \eqn{H_0=span\{1, sin(2\pi x),cos(2\pi x)\}};\cr
    \item linSinCos: \eqn{L=D^4+D^2}, \eqn{H_0=spac\{1, x, sin(x), cos(x)\}}.
}
}
\references{
Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.

Heckman, N and Ramsay, J. O. (2000). Penalised regression with model-based penalties.
To appear in Canadian Journal of Statisitcs.
}
\author{Chunlei Ke \email{chunlei\_ke@pstat.ucsb.edu} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{ssr}}
}
\examples{
\dontrun{
x<- seq(0,1, len=20)
lspline(x, type="exp", r=1.5)
}
}
\keyword{file}


