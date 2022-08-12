\name{Shrinkage}
\alias{shrink0}
\alias{shrink1}
\title{Calculate reproducing kernels for Stein shrinkage estimate}
\description{Return a matrix evaluating reproducing kernels for the discrete shrinkage towards zero or the mean estimate}
\usage{
shrink0(x, y=x)
shrink1(x, y=x)
}
\arguments{
\item{x}{a vector of numerical values or factor indicating different levels.}
\item{y}{a vector of numerical values or factor indicating different levels. Default is x.}
}
\value{
a matrix with the numbers of row and column equal to the length of x and y respectively. 
The \eqn{[i, j]} element is the reproducing kernel evaluated at the ith element of x and jth element of y.

\code{shink0} shrinks towards zero, and \code{shrink1} shinks towards the mean. 
}
\author{Chunlei Ke \email{chunlei_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{shrink0}},\code{\link{ssr}}
}
\examples{
\dontrun{
x<-rep(1:10,2)
shrink1(x)
}
}
\keyword{file}
