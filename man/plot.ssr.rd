\name{plot.ssr}
\alias{plot.ssr}
\title{Generate Diagnostic Plots for a ssr Object
}
\description{
Creates a set of plots suitable  for  assessing  a  fitted smoothing spline model of class \code{ssr}.
}
\usage{
\method{plot}{ssr}(x, ask=FALSE, ...)
}
\arguments{
   \item{x}{
    a \code{ssr} object.
   }
   \item{ask}{
   if TRUE, plot.ssr operates in interactive mode.
   }
   \item{...}{
   Other options used for plot, currently inactive.
  }
}


\details{
This function is a method for the  generic  function  plot for  class \code{ssr}.  
It can be invoked by calling plot for an object of the appropriate class, 
or  directly  by  calling plot.ssr regardless of the class of the object.

An appropriate x-y plot is produced to display  diagnostic plots.  These can be one or all of the following choices:
\itemize{
         \item Estimate of function with CIs 
         \item Residuals against Fitted values
         \item Response against Fitted values
         \item Normal QQplot of Residuals
}

The first plot of estimate of function with CIs is only useful for univariate smoothing spline fits.      

When ask=TRUE, rather than produce  each  plot  sequentially,  plot.ssr  displays a menu listing all the plots that can be produced. If the menu is not desired but a  pause  between plots  is  still  wanted  one  must  set par(ask=TRUE) before
invoking this command with argument ask=FALSE.
}
\author{Chunlei Ke \email{chunlei_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\seealso{
\code{plot}, \code{ssr}, \code{predict.ssr}
}
\examples{
\dontrun{library(MASS)}
\dontrun{fit1<- ssr(accel~times, data=mcycle, scale=TRUE, rk=cubic(times))}
\dontrun{plot(fit1,ask=TRUE)}
}
\keyword{file}
