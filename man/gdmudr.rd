\name{gdmudr}
\alias{gdmudr}
\title{Interface of dbmdr, dbimdr, dgmdr, dpmdr in GRKPACK.}
\description{
To calculate a spline estimate with multiple smoothing parameters for non-Gaussian data
}
\usage{
gdmudr(y, q, s, family, vmu = "v", varht = NULL, 
    init = 0, theta = NULL, tol1 = 0, tol2 = 0, prec1 = 1e-06, 
    maxit1 = 30, prec2 = 1e-06, maxit2 = 30) 
}
\arguments{
   \item{y}{
   a numerical vector representing the response, or a matrix of two columns for binomial data with the first column as the largest possible counts and the second column as the counts actually obsered.
    }
   \item{q}{
   a list, or an array, of square matrices of the same order as the length of y, which are the reproducing kernels evaluated at the design points.
   }
   \item{s}{
   the design matrix of the null space \eqn{H_0} of size (length-of-y,\eqn{dim(H_0)}), with elements equal to the bases of \eqn{H_0} evaluated at design points.
   }
   \item{family}{
   a string specifying the family of distribution. Families  supported  are  "binary", "binomial", "poisson" and "gamma" for Bernoulli, binomial, poisson, and gamma distributions respectively. Canonical links are used except for Gamma family where log link is used.
  }
   \item{vmu}{
   a character string specifying a method for choosing the smoothing  parameter.  "v", "m" and "u" represent GCV, GML and UBR respectively. 
"u\eqn{\sim}{~}", only used for non-Gaussian family, specifies UBR with estimated variance. Default is "v".
  }
  \item{varht}{
  needed only when vmu="u", which gives the fixed variance in calculation of the UBR function. Default is 1.0.
  }
  \item{init}{
   an integer of 0 or 1 indicating if initial values are provided for theta. If init=1, initial values are provided using theta. Default is 0.
  }
  \item{theta}{
  If `init=1', theta includes intial values for smoothing parameters. Default is NULL.
  }
  \item{tol1}{
  the tolerance for elements of w's. Default is 0.0 which sets to square of machine precision. 
  }
  \item{tol2}{
  tolerance for truncation used in `dsidr'. Default is 0.0 which sets to square of machine precision.
  }
  \item{prec1}{
  precision requested for the minimum score value, where precision is the weaker of the absolute and relative precisions. Default is 1e-06.
  }
  \item{maxit1}{
   maximum number of iterations allowed for DMUDR subroutine. Default is 30.
  }
  \item{prec2}{
  precision requested for stopping the iteration. Default is \eqn{1e-06}.
  }
  \item{maxit2}{
  maximum number of iterations allowed for the iteration in GRKPACK. Default is 30.
  }
}
\value{
  \item{info}{
  an integer that provides error message. info=-1 indicates dimension error, 
info=-2 idicates \eqn{F_{2}^{T} Q_{*}^{theta} F_{2} !>= 0}, info=-3 indicates tuning parameters are out of scope, info=-4 indicates dmudr fails to converge within maxit1 steps, info=-5 indicates dmudr fails to find a reasonable descent direction, info=-6 indicates GRKPACK fails to converge within maxit2 steps, info=-7 indicates there are some w's equals to zero, 
info>0 indicates the matrix S is rank deficient with \eqn{info=rank(S)+1}.
  }
  \item{fit}{
  estimate of the function at design points.
  }
  \item{c}{
  estimates of c.
  }
  \item{ d}{
  estimates of d.
  }
  \item{resi}{
   vector of working residuals.
  }
  \item{varht}{
  estimate of dispersion parameter.
  }
  \item{theta}{
  estimates of parameters \eqn{log10(theta)}. 
  }
  \itme{nlaht}{
  the estimate of \eqn{log10(nobs*lambda)}.
  }
  \item{score}{
  the minimum GCV/GML/UBR score at the estimated smoothing parameters. 
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
  \item{nq}{
  length(rk), number of reproducing kernels.
  }
  \item{s,q,y,init,maxit2}{
  changed from the inputs.
  }
}
\references{
Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.

Wang, Y. (1997). GRKPACK: Fitting Smoothing Spline ANOVA Models for Exponential Families. Communications in Statistics: Simulation and Computation, 24: 1037-1059.
}
\author{Chunlei Ke \email{chunlei\_ke@pstat.ucsb.edu} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{dsidr}}, \code{\link{dmudr}}, \code{\link{gdsidr}}, \code{\link{ssr}}
}
\keyword{file}
