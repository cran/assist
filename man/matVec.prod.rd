\name{matVec.prod}
\alias{matVec.prod}
\title{Product of a Matrix and a Vector}
\description{
 This function is used internally.
}
\usage{
matVec.prod(mat, vec, left = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{a matrix.}
  \item{vec}{a numeric vector of the same length as that of \code{mat}'s row or column}
  \item{left}{a logical value. If TRUE, multiply the matrix from the left. Otherwise it is from the right. Defaul is TRUE }
}
\details{
 If left=TRUE, perform \eqn{t(\mbox{vec})\%*\%\mbox{mat}}. Or else, perform
\eqn{\mbox{mat}\%*\%\mbox{vec}}.
}
\value{
 Returned is a vector.
}
\keyword{file}

