\name{deviance.ssr}
\alias{deviance.ssr}
\title{Model Deviance} 
\description{
Extract deviance from a fitted ssr object
}

\usage{
\method{deviance}{ssr}(object,residuals=FALSE, ...)
}
\arguments{
  \item{object}{a fitted \code{ssr} object}.
  \item{residuals}{
 a logical value. If 'TRUE', deviance residuals are returned. If 'FALSE', the sum of deviance residuals squares is returned. Default is FALSE.}
  \item{\dots}{other arguments, currently unused.}
}
\details{
       This is a method for the function \code{deviance}  for  objects
       inheriting from class \code{ssr}.
}
\author{Chunlei Ke \email{chunlei\_ke@pstat.ucsb.edu} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{ \code{\link{ssr}}}
\keyword{file}
