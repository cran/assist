\name{rkEval}
\alias{rkEval}
\title{Evaluate Reproducing Kernels}
\description{
This function evaluates a reproducing kernel or a list of reproducing kernels.
}
\usage{
rkEval(rk, g, g2)
}
\arguments{
  \item{rk}{a list of reproducing kernel expressions.}
  \item{g}{a data frame}
  \item{g2}{another data frame that may be different from \code{g}. Both data frames shall contain the variables appearing in \code{rk}.}
}
\details{
This function is used internally.
}
\value{
 a matrix or a list of matrices.
}
\keyword{file}


