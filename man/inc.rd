\name{inc}
\alias{inc}
\title{
Fit a Monotone Curve Using a Cubic Spline
}
\description{
Return a spline fit of a increasing curve.
}
\usage{
inc(y, x, spar = "v", limnla = c(-6, 0), grid = x, prec = 1e-06, maxit = 50, verbose = F)
}
\arguments{
  \item{y}{
a vecetor, used as the response data
}
  \item{x}{
a vector, used as the covariate. Assume an increasing relationshop of \code{y} on \code{x}
}
  \item{spar}{
a character string specifying a method for choosing the smoothing parameter. "v", "m" and "u" represent 
GCV, GML and UBR respectively. Default is "v" for GCV}
  \item{limnla}{
a vector of length one or two, specifying a search range for log10(n*lambda), where lambda is the smoothing 
parameter and n is the sample size. If it is a single value, the smoothing parameter will be fixed at this value. }
  \item{grid}{
a vector of \code{x} used to assess convergence. Default is \code{x}
}
  \item{prec}{
a numeric value used to assess convergence. Default is 1e-6
}
  \item{maxit}{
an integer represeenting the maximum iterations. Default is 50.
}
  \item{verbose}{
an optional logical value. If `TRUE', detailed iteration results are displayed. Default is "FALSE" 
}
}
\details{
This function is to fit a increasing fucntion to the data. The monotone function is expressed as integral of an unknown function that a cubic spline is used to estimate.
}
\value{
a split fit together with the convergence information
}
\author{Yuedong Wang \email{yuedong@pstat.ucsb.edu} and Chunlei Ke \email{chunlei\_ke@yahoo.com} }
\seealso{
\code{ssr}
}
\keyword{file }