\name{exponen}
\alias{exponen}
\title{Calculate Reproducing Kernels for L-spline based on Exponetial Functions}
\description{
Return a matrix evaluating a reproducing kernel for L-spline based on a differential 
operator involving exponential functions at observed points.
}
\usage{
exponen(s, t = s, r = 1)
}

\arguments{
  \item{s}{ a vector of values, at which the kernels are evaluated.  }
  \item{t}{an optional vector in [0, 1]. Default is the same as s.}
  \item{r}{a numeric value (integer) specifying the order. Default is 1.}
}
\details{
  An exponential L-spline is based on a linear differential operator 
\eqn{L=rD+D^2} with \eqn{H_0=span\{1,exp(-rx)\}} and \eqn{r>0}.
}
\value{
 a matrix with the numbers of row and column equal to the lengths
     of s and t respectively. The [i, j] element is the reproducing
     kernel evaluated at (s[i], t[j]).

}
\keyword{file}

