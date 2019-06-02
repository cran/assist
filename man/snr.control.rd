\name{snr.control}
\alias{snr.control}
\title{
Set Control Parameters for snr
}
\description{
Control parameters supplied in the function call replace 
the defaults to be used in calling \code{snr}.
}
\usage{
snr.control(rkpk.control = list(job = -1, tol = 0, init = 0, limnla = c(-10, 
    0), varht = NULL, theta = NULL, prec = 1e-06, maxit = 30), 
    nls.control = list(returnObject = TRUE, maxIter = 5), incDelta = 0.001, 
    prec.out = 0.001, maxit.out = 30, converg = "COEF", method = "GN", 
    backfit = 5) 
}
\arguments{
	\item{rkpk.control}{
	a optional list of control parameters for dsidr or dmudr to estimate the unknown
	functions. Default is "list(job = -1, tol = 0, init = 0, limnla = c(-10, 
    	0), varht = NULL, theta = NULL, prec = 1e-06, maxit = 30)".
	}
	\item{nls.control}{
	a list of control parameters for the nonlinear regression step, 
	the same as gnlsControl. Default is "list(returnObject = TRUE, maxIter = 5).
	}
        \item{incDelta}{the incremental value to be used to calculate derivatives for the unknown functions. Default is 0.001} 
	\item{prec.out}{
	tolerance for convergence criterion. Default is 0.0001.
	}
	\item{maxit.out}{
	maximum number of iterations for the algorithm. Default is 30.
	}
	\item{converg}{
	an optional character, with possible values \code{COEF} and \code{PRSS}, specifying the convergence 
	criterion to be used. \code{COEF} uses the change of estimate of parameters and functions to
	assess convergence, and \code{PRSS} uses penalized residual sums of squares. Default is \code{COEF}.
	}
	\item{method}{ an optional string of value either \code{GN} for Gauss-Newton or \code{NR} for Newton-Raphson 
           iteration methods to estimate the unknown functions. Default is \code{GN}.}
	\item{backfit}{ an integer to set the number of backfitting iterations inside the loop. Default is 5}.
}
\value{
returned is a list includes all re-seted control parameters.
}
\author{Chunlei Ke \email{chunlei\_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}.}
\seealso{
\code{\link{snr}}, \code{\link{dsidr}}, \code{\link{dmudr}}
}
\examples{
## use Newton-Raphson iteration and only a single backfitting
\dontrun{
snr.control(method="NR", backfit=1)
}
}
\keyword{file}
