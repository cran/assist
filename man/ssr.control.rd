\name{ssr.control}
\alias{ssr.control}
\title{Set Control Parameters for `ssr'}
\description{
The values supplied in the function call replace the defaults and a
list with all possible arguments is returned. The returned list is
used as the `control' argument to the `ssr' function.}
\usage{
ssr.control(job=-1, tol=0.0, init=0.0, theta, prec=1e-06, 
  maxit=30, tol.g=0.0, prec.g=1e-06, maxit.g=30)}
\arguments{
\item{job}{
an integer representing the optimization method used to find the smoothing parameter. 
The options are job=-1: golden-section search on (limnla(1), limnla(2)); 
job=0: golden-section search with interval specified automatically; 
job >0: regular grid search on \eqn{[limnla(1), limnla(2)]} with \#(grids) = job + 1. 
Default is -1. This is only applicable to smoothing spline model with 
a single smoothing parameter.
}
\item{tol}{
tolerance for truncation used in `dsidr' or `dmudr'. Default is 0.0 which sets to square of machine precision.
}
\item{init}{
init=0 means no initial values are provided for smoothing parameters theta; init=1 means initial values are provided for the theta. 
Default is 0. This option is only applicable to smoothing spline models with multiple smoothing parameters.
}
\item{theta}{
If init=1, theta includes intial values for smoothing parameters. Default is NULL. This is only applicable to smoothing spline models with multiple smoothing parameters.
}
\item{prec}{
precision requested for the minimum score value in `dmudr', where precision is the weaker of the absolute and relative precisions. 
Default is 1e-06. This is only applicable to smoothing spline models with multiple smoothing parameters.
}
\item{maxit}{
maximum number of iterations allowed in `dmudr'. Default is 30. This is only applicable to smoothing spline model with multiple smoothing parameters.
}
\item{tol.g}{
the tolerance for elements of w's in GRKPK. Default is 0.0 which means using the machine precision. This is only applicable to generalized spline smoothing.
}
\item{prec.g}{
precision for stopping the iteration in GRKPK. Default is 1e-06. This is only applicale to generalized spline smoothing.
}
\item{maxit.g}{
maximum number of iterations allowed for the iteration in GRKPACK. Default is 30. This is only applicale to generalized spline smoothing.}
} 
\value{
a list with components for each of the possible arguments.
}
\seealso{
\code{\link{ssr}}
}
\example{
# use regular grid seach method with 100 grid points
ssr.control(job=99)
}
\keyword{file}

