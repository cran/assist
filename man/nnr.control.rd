\name{nnr.control}
\alias{nnr.control}
\title{Set Control Parameters for nnr}
\description{
Control parameters supplied in the function call replace 
the defaults to be used in calling \code{nnr}.
}
\usage{
nnr.control(job = -1, tol = 0, max.iter = 50, init = 0, limnla = c(-10, 
    0), varht = NULL, theta = NULL, prec = 1e-06, maxit = 30, 
    method = "NR", increment = 1e-04, backfit = 5, converg = "coef", 
    toler = 0.001)   
}
\arguments{
\item{job}{
  an integer representing the optimization method used to find the smoothing parameter. 
The options are job=-1: golden-section search on (limnla(1), limnla(2)); 
job=0: golden-section search with interval specified automatically; 
job >0: regular grid search on  [limnla(1), limnla(2)] with \#(grids) = job + 1. 
Default is -1. 
  }
\item{ tol }{
  tolerance for truncation used in `dsidr'. Default is 0.0, which sets to square of machine precision.
  }
\item{max.iter}{maximum number of iterations allowed for the Gauss-Newton/Newton-Raphson iteration.}
\item{init}{
  an integer of 0 or 1 indicating if initial values are provided for theta. If init=1, initial values are provided using theta. Default is 0.
  }
\item{limnla}{
  a vector of length 2, specifying a search range for the  n times smoothing parameter on log10 scale. Default is (-10, 0).
  }
 \item{ varht}{
  needed only when vmu="u", which gives the fixed variance in calculation of the UBR function. Default is NULL.
  }
 \item{theta}{
  If `init=1', theta includes intial values for smoothing parameters. Default is NULL.
  }
 \item{ prec}{
  precision requested for the minimum score value, where precision is the weaker of the absolute and relative precisions. Default is 1e-06.
  }
  \item{ maxit}{
   maximum number of iterations allowed. Default is 30.
   }
\item{method}{a character string specifying a method for iterations, "GN" for Gauss-Newton and "NR" for Newton-Raphson. Default is "GN".}
\item{increment} {specifies a small value as increment to calcuate derivatives. Default is 1e-04.}
\item{backfit}{ an integer representing the number of backfitting iterations for multiple functions. Default is 5.}
\item{converg}{	an optional character, with possible values "coef" and "ortho", specifying the convergence 
	criterion to be used. "coef" uses the change of estimate of parameters and functions to
	assess convergence, and "ortho" uses a criterion similar to the relative offset used in nls. Default is "coef".
	}
\item{toler}{tolerance for convergence of the algorithm. Default is 0.001.}
}
\value{returned is a list includes all re-seted control parameters.
}
\author{Chunlei Ke \email{chunlei\_ke@pstat.ucsb.edu} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{nnr}}, \code{\link{dsidr}},\code{\link{dmudr}}
}
\examples{
## use Newton-Raphson 
nnr.control(method="NR")
}
\keyword{file}
