\name{ident}
\alias{ident}
\title{Scaling a Vector}
\description{
Perform standarization of vector relative to another.
}
\usage{
ident(x, y = x)
}
\arguments{
  \item{x}{a numeric vector, matrix or data frame}
  \item{y}{an optional numeric vector, matrix or data frame. Default is x.}
}
\details{
Scale \code{y} based on \code{x} component by component. For example, if both are a matrix, 
\eqn{y[,i]=(y[,]-min(x[,i]))/(max(x[,i])-min(x[,i]))}.
}
\value{
 a scaled \code{y}.
}
\keyword{file}


