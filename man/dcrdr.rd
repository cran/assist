\name{dcrdr}
\alias{dcrdr}

\title{Interface to Fortran Subroutine dcrdr}
\description{
  Calculate some matrix operations needed to construct Bayesian confidence intervals
}
\usage{
dcrdr(rkpk.obj, r)
}
\arguments{
  \item{rkpk.obj}{an object returned from calling dsidr}
  \item{r}{a matrix to evaluate reproducing kernels on grid points}
}
\value{
See the document for the corresponding Fortran subroutine.
}
\keyword{file}


