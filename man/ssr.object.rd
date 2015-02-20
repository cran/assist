\name{ssr.object}
\alias{ssr.object}
\title{ A fitted ssr Object}
\description{
An object returned by the \code{ssr} function, inheriting from class \code{ssr},
and representing a fitted smoothing spline model. Objects of this
class have methods for the generic functions \code{predict}, \code{print} and
\code{summary}. }
\value{
The following components must be included in a legitimate \code{ssr} object: 
\item{call}{a list containing an image of the \code{ssr} call that produced the object}
\item{coef}{estimated coefficients for the spline estimate}
\item{lambda}{a vector representing the estimate smoothing parameters}
\item{fitted}{fitted values of the unknown mean function}
\item{family}{the distribution family used}
\item{cor.est}{estiamted parameters, if any, in corMatrix}
\item{var.est}{estiamted parameters, if any, in varFunc}
\item{s}{design matrix extracted from \code{formula}}
\item{q}{a list of matrices representing reproducing kernels evaluated at design points.}
\item{residuals}{working residuals from the fit. }
\item{df}{equivalent degrees of freedom. It is calculated as the trace of the hat matrix.}
\item{weight}{a matrix representing the covariance matrix. It is NULL for iid data.}
\item{rkpk.obj}{an object representing fits from dsidr/dmudr/gdsidr/gdmudr. See help files 
        for dsidr/dmudr/gdsidr/gdmudr for more details.}
\item{scale}{a logical value, specifying if scaling is used.}
}

\author{Chunlei Ke \email{chunlei\_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}

\seealso{
\code{\link{ssr}}, \code{\link{predict.ssr}}, \code{\link{summary.ssr}}, 
\code{\link{plot.ssr}}, \code{\link{dsidr}}, \code{\link{dmudr}}, \code{\link{gdsidr}},
\code{\link{gdmudr}}
}
\keyword{file}
