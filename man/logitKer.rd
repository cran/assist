\name{logitKer}
\alias{logitKer}
\title{Calculate Reproducing Kernels for a L-spline Based on the Logit Function }
\description{
Return a matrix evaluating a reproducing kernel for L-spline based on a differential 
operator with the logistic function in the null space.
}
\usage{
logitKer(s, t = s)
}
\arguments{
\item{s}{a vector of values, at which the kernels are evaluated.  }
\item{t}{an optional vector. Default is the same as s.}
}
\details{
 An exponential L-spline is based on a linear differential operator 
  \eqn{L=D-1/(1+e^t)} with \eqn{H_0=span\{e^t/(1+e^t)\}}.  
}
\value{
 a matrix with the numbers of row and column equal to the lengths
     of s and t respectively. The [i, j] element is the reproducing
     kernel evaluated at (s[i], t[j]).
}
\keyword{file}

