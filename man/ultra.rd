\name{ultrasound}
\docType{data}
\alias{ultrasound}
\title{
Ultrasound imaging of the tongue shape
}
\description{
The 'ultrasound' data frame has 1,215 rows and 4 columns of data from an ultrasound experiment
}
\usage{
data(ultrasound)
}
\section{Format}{
The data frame contains the following columns:

height a numeric vector of toungue height in mm

length a numeric vector of toungue length in mm

time a numeric vector of time in ms

env a factor with three levels: 1 2 and 3 for environment '2words', 'cluster' and 'Schwa' respectively
}
\details{
A Russian speaker produced the consonant sequence, /gd/, in three different linguistic environments: '2words', 'cluster' 
and 'Schwa', with three replications for each environment. 15 points from each of 9 slices of toungue curves separated by 30 ms (milliseconds) 
are extracted. Therefore, in total there are 15*9*3*3=1,215 observations.
}
\section{Source}{
Phonetics-Phonology Lab of New York University. 
}
\section{references}{
Davidson, L. (2006). Comparing Tongue Shapes from Ultrasound Imaging Using Smoothing Spline Analysis of Variance. Journal of the Acoustical Society of America 120, 407-415. 
}
\keyword{datasets}