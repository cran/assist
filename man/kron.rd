\name{kron}
\alias{kron}
\title{Calculate reproducing kernels for one-dimensional space}
\description{Return a matrix evaluating reproducing kernels for the one-dimensional space usually spanned by a vector}
\usage{ kron(x,y=x)
}
\arguments{
   \item{x}{a vector or a list of numerical values which spans the one-dimensional space.}
   \item{y}{a vector or a list of numerical values. Default is x.}
}
\value{
a matrix with the numbers of row and column equal to the length of x and y respectively. 
The [i, j] element is the reproducing kernel evaluated at the ith element of x and jth element of y.
}
\author{Chunlei Ke \email{chunlei_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{\link{kronecker}},\code{\link{ssr}}
}
\examples{
\dontrun{
x<-runif(10)
kron(x)
}
}
\keyword{file}
