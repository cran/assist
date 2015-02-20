\name{bond} 
\docType{data}
\alias{bond}            
\title{Treasury and GE bonds}
\description{The \code{bond} data frame has 1234 rows and 5 columns of data derived from 144 General Electronic Company bonds and 78 Treasury bonds.}
\usage{
data(bond)
}     
\section{Format}{
The data frame contains the following columns:

name a vector of index for individual bond

price a numeric vector of current price

time a numeric vector of future time points at which the payments are made

payment a numeric vector of future payments
 
type a vector of character strings identifying the groups, "govt" or "ge", which the individual bonds belong to.
}
\section{Source}{
Bloomberg
}
\section{references}{
Chunlei Ke and Yuedong Wang (2004), Nonlinear Nonparametric Regression Models. Journal of the American Statistical Association 99, 1166-1175.
}
\keyword{datasets}



