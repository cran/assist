\name{xyplot2}
\alias{xyplot2}
\title{Extension of XYPLOT}
\description{
Extend \code{xyplot} to superpose one or more symbols to each panel.
}
\usage{
xyplot2(formula, data, type = "l", ...)
}
\arguments{
  \item{formula}{a two-sided formula as accepted in \code{xyplot}}
  \item{data}{a list of data frames. Each component shall be able to evaluate the vatiables appearing in \code{formula}}
  \item{type}{a vector of characters to indicate what type of plots are to draw. Default is line.}
  \item{\dots}{any options as accepted in \code{xyplot}}
}
\value{
On each panel, several plot types, the length of \code{data}, are superposed.
}
\keyword{file}


