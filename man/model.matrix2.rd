\name{model.matrix2}
\alias{model.matrix2}
\title{Extension of the Model Matrix Function}
\description{
This function is employed internally.
}
\usage{
model.matrix2(formula, data = list(), ...)
}
\arguments{
  \item{formula}{a formula as accepted in \code{model.matrix}.}
  \item{data}{a data frame. If missing, the data are to get from the search list. }
  \item{\dots}{other options as accepted in \code{model.matrix}}
}
\details{
  The only extension is when \code{formula} is \eqn{\sim}{~} 1. The returned matrix is a constant vector of the same length as the number of observations in \code{data}.
}
\value{
 A matrix, represent the model specified in \code{formula}.
}
\keyword{file}


