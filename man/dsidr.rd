\name{ dsidr }
\alias{dsidr}
\title{Interface of dsidr subroutines in RKPACK}
\description{
To calculate a spline estimate with a single smoothing parameter
}
\usage{dsidr(y, q, s=NULL, weight=NULL, vmu="v", varht=NULL, 
limnla=c(-10, 3), job=-1, tol=0)
}
\arguments{
  \item{y}{
   a numerical vector representing the response.
  }
  \item{q}{
   a square matrix of the same order as the length of y, with elements equal to the reproducing kernel evaluated at the design points.
  }
  \item{s}{
  the design matrix of the null space \eqn{H_0} of size (length(y),dim(\eqn{H_0})), 
with elements equal to the bases of \eqn{H_0} evaluated at design points. Default is NULL, representing an empty NULL space.}
  \item{weight}{
  A weight matrix for penalized weighted least-square: \eqn{(y-f)'W(y-f)+n\lambda J(f)}. Default is NULL for iid random errors.
  }
 \item{vmu}{
  a character string specifying a method for choosing the smoothing  parameter.  "v", "m" and "u" represent GCV, GML and UBR respectively. 
"u\eqn{\sim}{~}", only used for non-Gaussian family, specifies UBR with estimated variance. Default is "v".
 }
 \item{varht}{
  needed only when vmu="u", which gives the fixed variance in calculation of the UBR function. Default is NULL.
  }
  \item{limnla}{
  a vector of length 2, specifying a search range for the  n times smoothing parameter on \eqn{log10} scale. Default is \eqn{(-10, 3)}.
  }
  \item{job}{
  an integer representing the optimization method used to find the smoothing parameter. 
The  options are job=-1: golden-section search on (limnla(1), limnla(2)); 
job=0: golden-section search with interval specified automatically; 
job >0: regular grid search on  \eqn{[limnla(1), limnla(2)]} with the number of grids = job + 1. Default is -1. 
  }
  \item{tol }{
  tolerance for truncation used in `dsidr'. Default is 0.0, which sets to square of machine precision.
  }
}
\value{
  \item{info}{
   an integer that provides error message. info=0 indicates normal termination, info=-1 indicates dimension error, info=-2 indicates
 \eqn{F_{2}^{T} Q F_{2} !>= 0}, info=-3 indicates vmu is out of scope, and info>0 indicates the matrix S is rank 
deficient with info=rank(S)+1. 
   }
  \item{ fit}{
  fitted values.
  }
  \item{c}{
  estimates of c.
  }
  \item{d}{
  estimates of d.
  }
  \item{resi}{
  vector of residuals.
  }
  \item{varht}{
  estimate of variance.
  }
  \item{nlaht}{
  the estimate of log10(nobs*lambda).
  }
  \item{limnla}{
  searching range for nlaht. 
  }
  \item{score}{
  the minimum GCV/GML/UBR score at the estimated smoothing parameter. When job>0, it gives a vector of GCV/GML/UBR functions evaluated at regular grid points.
  }
  \item{df}{
  equavilent degree of freedom.
  }
  \item{nobs}{
  length(y), number of observations.
  }
  \item{nnull}{
  dim(\eqn{H_0}), number of bases.
  }
  \item{s,qraux,jpvt}{
   QR decomposition of S=FR, as from Linpack `dqrdc'.
   }
  \item{q}{
   first dim(\eqn{H_0}) columns gives \eqn{F^{T} Q F_{1}}, and its bottom-right corner gives tridiagonalization of \eqn{F_{2}^{T} Q F_{2}}.
   }
}
\references{
Gu, C. (1989). RKPACK and its applications: Fitting smoothing spline models. Proceedings of the Statistical Computing Section, ASA, 42-51.

Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.
}
\author{Chunlei Ke \email{chunlei_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{dmudr}}, \code{\link{gdsidr}}, \code{\link{gdmudr}}, \code{\link{ssr}}
}
\keyword{file}
