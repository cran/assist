\name{TXtemp}
\docType{data}
\alias{TXtemp}
\title{Texas Historical Climate Data}
\description{The data frame \code{TXtemp}, obtained from the Carbon Dioxide Information and Analysis Center at Oak Ridge National Laboratory, has 17280 rows and 6 columns of data representing monthly temperature records for stations in Texas.}

\usage{
data(TXtemp)
}
\format{
The data frame contains the following columns:

stacode a numeric vector of the unique station code formed by combining the two-digit state number 
[state numbers range from 1 to 48] and the four-digit station number (values range from 0008 to 9933);                   

lat, long numeric vectors identifying the lattitudes and longitudes of the stations in decimal degree.

year a numeric vector comprising the year for the records

month a numeric vector of values 1 to 12, represeting the month for the data

mmtemp a numeric vector of monthly average temperature in Fahrenheit scale.
}

\details{The data set was extracted from a large national historical climate data, containing data for 48 stations in Texas from 1961 to 1990. Monthly temperature records as well as the latitude and longitude for each station were available.

Of note, the missing values were coded as -99.99.
}

\source{
Data are downloadable from \url{http://cdiac.ornl.gov/ftp/ndp019/}
}
\keyword{datasets}
