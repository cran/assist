\name{gdsidr}
\alias{gdsidr}
\title{
Interface of dbsdr, dbisdr, dgsdr, dpsdr in GRKPACK.
}
\description{
To calculate a spline estimate with single smoothing parameter for non-Gaussian data.
}
\usage{
gdsidr(y, q, s, family, vmu="v", varht=NULL, limnla=c(-10, 3), 
maxit=30, job=-1, tol1=0, tol2=0, prec=1e-06)
}
\arguments{
 \item{ y}{
  a numerical vector representing the response, or a matrix of two columns for binomial data with the first column as the largest possible counts and the second column as the counts actually obsered.
  }
 \item{q}{
  a square matrix of the same order as the length of y, with elements equal to the reproducing kernel evaluated at the design points.
  }
  \item{s}{
  the design matrix of the null space \eqn{H_0} of size (length-of-y,dim(\eqn{H_0})), with elements equal to the bases of \eqn{H_0} evaluated at design points.
  }
  \item{family}{
  a string specifying the family of distribution. Families  supported  are  "binary", "binomial", "poisson" and "gamma" for Bernoulli, binomial, poisson, and gamma distributions respectively. Canonical links are used except for Gamma family where a log link is used.
  }
  \item{vmu}{
  a character string specifying a method for choosing the smoothing  parameter.  "v", "m" and "u" represent GCV, GML and UBR respectively. "u\eqn{\sim}{~}", only used for non-Gaussian family, specifies UBR with estimated variance. Default is "v".
  }
  \item{varht}{
  needed only when vmu="u", which gives the fixed variance in calculation of the UBR function. Default is 1.0.
  }
  \item{limnla}{
  a vector of length 2, specifying a search range for the  n times smoothing parameter on log10 scale. Default is (-10, 3).
  }
  \item{maxit}{
  maximum number of iterations allowed for the iteration in GRKPACK.
  }
  \item{job}{
  an integer representing the optimization method used to find the smoothing parameter. 
The  options are job=-1: golden-section search on (limnla(1), limnla(2)); 
job=0: golden-section search with interval specified automatically; 
job >0: regular grid search on  [limnla(1), limnla(2)] with the number of grids = job + 1. Default is -1. 
  }
  \item{tol1}{
  the tolerance for elements of w's. Default is 0.0 which sets to square of machine precision. 
  }
  \item{tol2}{
  tolerance for truncation used in `dsidr'. Default is 0.0 which sets to square of machine precision.
  }
  \item{prec}{
  precision requested for stopping the iteration. Default is \eqn{1e-06}.
  }
}
\value{
  \item{info}{
   an integer that provides error message. info=0 indicates normal termination, info=-1 indicates dimension error, 
info=-2 indicates \eqn{F_{2}^{T} Q F_{2} !>= 0}, info=-3 indicates vmu is out of scope, info=-4 indicates the algorithm fails to converge at the maxiter steps, info=-5 indicates there are some w's equals to zero, and info>0 indicates the matrix S is rank deficient with info=rank(S)+1. 
   }
  \item{fit}{
   estimate of the function at design points.
   }
  \item{c}{
   estimates of c.
  }
  \item{d}{
  estimates of d.
  }
  \item{resi}{
   vector of working residuals.
   }
  \item{varht}{
  estimate of dispersion parameter.
  }
  \item{nlaht}{
  the estimate of \eqn{log10(nobs*lambda)}.
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
  length-of-y, number of observations.
  }
  \item{nnull}{
   \eqn{dim(H_0)}, number of bases.
  }
  \item{s,qraux,jpvt}{
  QR decomposition of S=FR, as from Linpack `dqrdc'.
  }
  \item{q}{
   first \eqn{dim(H_0)} columns gives \eqn{F^{T} Q F_{1}}, and its bottom-right corner gives tridiagonalization of \eqn{F_{2}^{T} Q F_{2}}.
   }
}
\references{
Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.

Wang, Y. (1997). GRKPACK: Fitting Smoothing Spline ANOVA Models for Exponential Families. Communications in Statistics: Simulation and Computation, 24: 1037-1059.
}
\author{Chunlei Ke \email{chunlei_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{dsidr}}, \code{\link{dmudr}}, \code{\link{gdmudr}}, \code{\link{ssr}}
}
\keyword{file}
