\name{diagComp}
\alias{diagComp}
\title{Form a Block Diagonal Matrix}
\description{
  This function is used internally
}
\usage{
diagComp(mat, vec)
}

\arguments{
  \item{mat}{a matrix of any order}
  \item{vec}{a vector of index indicating the number of columns for each block}
}
\value{
 Returned is a block diagnoal matrix. The number of blocks is equal to the length of \code{vec},
and each block has the same number of rows as \code{mat}. The blocks' columns are defined by \code{vec}.
}
\keyword{file}
