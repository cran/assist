\name{linSinCos}
\alias{linSinCos}
\title{Calculate Reproducing Kernels for linear periodic L-spline}
\description{
 Return a matrix evaluating a reproducing kernel for L-spline based on a differential 
operator with the linear and sine and cosine functions in the null space.
}
\usage{
linSinCos(s, t = s)
}

\arguments{
\item{s}{ a vector of values, at which the kernels are evaluated.  }
\item{t}{an optional vector. Default is the same as s.}
}
\details{
 An exponential L-spline is based on a linear differential operator 
 \eqn{L=D^4+D^2} with \eqn{H_0=spac\{1, x, sin(x), cos(x)\}}.  
}
\value{
 a matrix with the numbers of row and column equal to the lengths
     of s and t respectively. The [i, j] element is the reproducing
     kernel evaluated at (s[i], t[j]).
}
\keyword{file}

