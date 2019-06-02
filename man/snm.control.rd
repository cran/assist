\name{snm.control}
\alias{snm.control}
\title{
Set Control Parameters for snm
}
\description{
Control parameters supplied in the function call replace 
the defaults to be used in calling \code{snm}.
}
\usage{
snm.control(rkpk.control, nlme.control, prec.out=0.0005, 
  maxit.out=30, converg="COEF", incDelta)
}
\arguments{
\item{rkpk.control}{
a optional list of control parameters for dsidr or dmudr to estimate the unknown
functions.
}
\item{nlme.control}{
a list of control parameters for the nonlinear regression step, 
the same as nlmeControl. Default is \code{list(returnObject = T, maxIter = 5)}.
}
\item{prec.out}{
tolerance for convergence criterion. Default is 0.0005.
}
\item{maxit.out}{
maximum number of iterations for the algorithm. Default is 30.
}
\item{converg}{
an optional character, with possible values "COEF" and "PRSS", specifying the convergence 
criterion to be used. "COEF" uses the change of estimate of parameters and functions to
assess convergence, and "PRSS" uses penalized residual sums of squares. Default is "COEF".
}
\item{incDelta}{specifies a small value as increment to calcuate derivatives. Default is 0.001.
}
}
\value{
Returned is a list includes all re-seted control parameters.
}
\author{Chunlei Ke \email{chunlei\_ke@pstat.ucsb.edu} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{snm}}, \code{\link{dsidr}}, \code{\link{dmudr}}
}
\examples{
\dontrun{
## set maximum iteration to be 50
snm.control(maxit.out=50)
}
}
\keyword{file}

