\name{matDiag.prod}
\alias{matDiag.prod}
\title{Product of a Matrix and a Diagonal Matric formed from a Vector}
\description{
This is a function used internally.
}
\usage{
matDiag.prod(mat, vec, left = TRUE)
}
\arguments{
  \item{mat}{a matrix.}
  \item{vec}{a numeric vector of the same length as that of \code{mat}'s row or column}
  \item{left}{a logical value. If TRUE, multiply the matrix from the left. Otherwise it is from the right. Defaul is TRUE }
}
\details{
If left=TRUE, perform \eqn{\mbox{diag(vec)}\%*\%\mbox{mat}}. Or else, perform
\eqn{\mbox{mat}\%*\%\mbox{diag(vec)}}.
}
\value{
 a matrix is returned.
}
\keyword{file}


