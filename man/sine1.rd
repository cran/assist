\name{sine1}
\alias{sine1}
\title{ Calculate Reproducing Kernels for a Periodic L-spline without the Constant Term}
\description{
  Return a matrix evaluating a reproducing kernel for L-spline based on a differential 
operator with the constant, sine and cosine functions in the null space. 
}
\usage{
sine1(s, t = s)
}
\arguments{
\item{s}{ a vector of values, at which the kernels are evaluated.}
\item{t}{an optional vector. Default is the same as s.}
}
\details{
 An exponential L-spline is based on a linear differential operator 
  \eqn{L=D^2+(2*pi)^2} with \eqn{H_0=span\{1, sin(2*pi*x),cos(2*pi*x)\}}.  
}
\value{
 a matrix with the numbers of row and column equal to the lengths
     of s and t respectively. The [i, j] element is the reproducing
     kernel evaluated at (s[i], t[j]).
}
\keyword{file}





