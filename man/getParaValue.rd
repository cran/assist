\name{getParaValue}
\alias{getParaValue}
\title{Evaluate Parameters}
\description{
Get parameters' values based on data and estimates
}
\usage{
getParaValue(length.fix, length.random, matrix.fix, matrix.random, start.fixed, start.random, para, para.random, nobs, para.fixed, nobs.in)
}
\arguments{
  \item{length.fix}{a vector}
  \item{length.random}{a vector }
  \item{matrix.fix}{the model matrix for fixed effects }
  \item{matrix.random}{the model matrix for random effects}
  \item{start.fixed}{start/estimate values for fixed effects }
  \item{start.random}{start/extimate values for random effects}
  \item{para}{a list of all parameter symbols}
  \item{para.random}{a list of random effects symbols }
  \item{nobs}{an integer}
  \item{para.fixed}{a list of fixed effects symbols}
  \item{nobs.in}{lengths in each group}
}
\details{
This function is used internally.
}
\keyword{file}

