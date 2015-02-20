\name{horm.cort}
\docType{data}
\alias{horm.cort}
\title{Hormone Measurements of Cortisol}
\description{The \code{horm.cort} data frame has 425 rows and 4 columns of data representing measurement of cortisol on 36 individuals.}

\usage{
data(horm.cort)
}
\format{
The data frame contains the following columns:

ID a vector of index indicating individuals on whom measures were made.

time a numeric vector of time points of every 2 hours in 24 hours. The time is scaled into [0, 1].

type a vector of character strings identifying the groups, "normal", "depressed", or "cushing", which the individuals belong to.

conc cortisol concentration measurements in \eqn{log10} scale.
}

\details{Blood samples were collected every 2 hours for 24 hours from three group of healthy normal volunteers and volunteers with depresession and suchsing syndrome. They were analyzed for parameters that measure hormones of the hypothalamic-pituitary axix. Human circadian thythm is one of the research objective. In this data set, only measurements of cortisol concetration were included.}

\source{
This data set was extracted from a stress study conducted in the medical center of the University of Michigan.
}
\references{
Wang, Y. and Brown, M. B. (1996). A Flexible Model for Human Circadian Rhythms. Biometrics 52, 588-596.

Yuedong Wang, Chunlei Ke and Morton B. Brown (2003), Shape Invariant Modelling of Circadian Rhythms with Random Effects and Smoothing Spline ANOVA Decompositions. Biometrics, 59:804-812.
}
\keyword{datasets}
