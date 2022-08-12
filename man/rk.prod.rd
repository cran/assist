\name{rk.prod}
\alias{rk.prod}
\title{
Calculate product of reproducing kernels
}
\description{
Return a matix as product of reproducing kernels
}
\usage{
rk.prod(x, \dots)
}
\arguments{
	\item{x}{
	a matrix evaluating a reproducing kernel, or a vector.
	}
	\item{\dots}{
	optional lists of matrices evaluating reproducing kernels or vectors. All matrics
	must have the same dimensions. All vectors must have the same length. The length of
	each vector must equal to the column  and row numbers of each matrix.
	}
}
\value{
a matrix as the product of reproducing kernels. If one argument is a vector, a \code{kron}
kernel is constructed first.
}
\details{
The product of reproducing kernels is agian a reproducing kernel. In SS ANOVA, product
of reproduing kernels is often used to model interaction spline terms.
}
\references{
Gu, C. and Wahba, G. (1993a). Smoothing Spline ANOVA with component-wise Bayesian confidence intervals.
Journal of Computational and Graphical Statistics 55, 353--368.

Gu, C. and Wahba, G. (1993b). Semiparametric analysis of variance with tensor product thin plate
splines. JRSS B 55, 353--368. 
}
\author{Chunlei Ke \email{chunlei_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{\code{\link{kron}}, \code{\link{ssr}}}
\examples{
\dontrun{
x1<- 1:10/10
x2<- runif(10)
rk.prod(cubic(x1), periodic(x2))
}
}
\keyword{file}


