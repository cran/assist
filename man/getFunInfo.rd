\name{getFunInfo}
\alias{getFunInfo}
\title{Extract the Unknown Spline Function Sructures}
\description{
Used internally, this function gets the null spaces and reporducing kernels for the unknown functions
}
\usage{
getFunInfo(object)
}

\arguments{
  \item{object}{a list of spline structures. }
}
\details{
The spline structure for an unknown function consists of "nb", a one side formula specifying the null space, and "rk", an expression representing the repoducing kernel. \code{getFunInfo} breaks down the structures into components including names, arguments, bases of the null space and rk for each function.
}
\value{
A list of length 4 consists of a list of function names, basis formulae for the null spaces, reproducing kernels and a list of symbols of arguments.
}
\keyword{file}


