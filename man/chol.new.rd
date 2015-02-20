\name{chol.new}
\alias{chol.new}
\title{A Modified Cholesky Decomposition}
\description{
 Returned a matrix forming Cholesky Decomposition
}
\usage{
chol.new(Q)
}
\arguments{
  \item{Q}{a symmetric matrix, maybe non-positive.}
}
\details{
 This is used internally as an extension of \code{chol} that works on a positive matrix.
}
\value{
 A mtrix M suth that \eqn{XX^T=Q}.
}
\seealso{\code{\link{chol}}}
\keyword{file}


