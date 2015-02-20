\name{dmudr}
\alias{dmudr}
\title{
Interface of dmudr subroutine in RKPACK
}
\description{
To calculate a spline estimate with multiple smoothing parameters
}
\usage{
dmudr(y, q, s, weight = NULL, vmu = "v", theta = NULL, varht = NULL, 
    tol = 0, init = 0, prec = 1e-06, maxit = 30)
}
\arguments{
 \item{y}{
  a numerical vector representing the response.
  }
  \item{q}{
  a list, or an array, of square matrices of the same order as the length of y, which are the reproducing kernels evaluated at the design points.
  }
  \item{s}{
  the design matrix of the null space \eqn{H_0} of size (length-of-y,\eqn{dim(H_0)}), with elements equal to the bases of \eqn{H_0} evaluated at design points.
  }
  \item{weight}{
   a weight matrix for penalized weighted least-square: \eqn{(y-f)'W(y-f)+n\lambda J(f)}. Default is NULL for iid random errors.
  }
  \item{vmu}{
  a character string specifying a method for choosing the smoothing  parameter.  "v", "m" and "u" represent GCV, GML and UBR respectively. "u\eqn{\sim}{~}", only used for non-Gaussian family, specifies UBR with estimated variance. Default is "v".
  }
   \item{theta}{
  If `init=1', theta includes intial values for smoothing parameters. Default is NULL.
  }
  \item{varht}{
  needed only when vmu="u", which gives the fixed variance in calculation of the UBR function. Default is NULL.
  }
  \item{tol}{
  the tolerance for truncation in the tridiagonalization. Default is 0.0.
  }
  \item{init}{
  an integer of 0 or 1 indicating if initial values are provided for theta. If init=1, initial values are provided using theta. Default is 0.
  }
  \item{prec}{
  precision requested for the minimum score value, where precision is the weaker of the absolute and relative precisions. Default is \eqn{1e-06}.
  }
  \item{maxit}{
   maximum number of iterations allowed. Default is 30.
   }
}
\value{
\item{info}{
  an integer that provides error message. info=-1 indicates dimension error, 
info=-2 indicates \eqn{F_{2}^{T} Q_{*}^{\theta} F_{2} !>= 0}, info=-3 indicates tuning parameters are out of scope, info=-4 indicates fails to converge within maxite steps, info=-5 indicates fails to find a reasonable descent direction, info>0 indicates the matrix S is rank deficient with \eqn{info=rank(S)+1}.
  }
   \item{fit}{
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
  \item{theta}{
   estimates of parameters \eqn{log10(\theta)}. }
  \item{nlaht}{
  the estimate of \eqn{log10(nobs*\lambda)}.
  }
  \item{score}{
  the minimum GCV/GML/UBR score at the estimated smoothing parameters. 
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
  \item{nq}{
  length(rk), number of reproducing kernels.
  }
  \item{s,q,y}{
  changed from the inputs.
  }
}
\references{
Gu, C. (1989). RKPACK and its applications: Fitting smoothing spline models. Proceedings of the Statistical Computing Section, ASA, 42-51.

Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59
}
\author{Chunlei Ke \email{chunlei\_ke@pstat.ucsb.edu} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{dsidr}}, \code{\link{gdsidr}}, \code{\link{gdmudr}}, \code{\link{ssr}}
}
\keyword{file}
