\name{bdiag}
\alias{bdiag}
\title{Construct a Block Diagonal Matrix}
\description{
  Return a block diagonal matrix formed from the input list of matrices
}
\usage{
bdiag(x)
}
\arguments{
  \item{x}{a list of matrices}
}
\value{
 Returned is a matrix of the form diag(x1, \dots, xn) where n is the length of the list.
}
\keyword{file}
