\name{seizure}
\docType{data}
\alias{seizure}            
\title{IEEG segments from a seizure patient}
\description{
The 'seizure' data frame has 60,000 rows and 3 columns of data from an IEEG time series
}
\usage{
data(seizure)
}     
\section{Format}{
The data frame contains the following columns:

t a numeric vector of the observation number

base a numeric vector of the baseline segment

preseizure a numeric vector of the segment right before a seizure
}
\details{
The baseline segment contains 5-minute IEEG signal
extracted at least four hours before the seizure's onset.
The preseizure segment contains 5-minute IEEG signal
right before a seizure's clinical onset. The sampling rate 
of the IEEG signal is 200 observations per second.
Therefore there are 60,000 time points in each segment.
}
\section{Source}{
D'Alessandro, M., Vachtsevanos, G., Esteller, R., Echauz, J. 
and Litt, B. (2001). A Generic Approach to Selecting the 
Optimal Feature for Epileptic Seizure Prediction.
IEEE International Meeting of the Engineering in Medicine and Biology Society.
}
\section{references}{
Qin, L. and Wang, Y. (2008), Nonparametric Spectral Analysis
With Applications to Seizure Characterization Using EEG Time Series. Annals of Applied Statistics 2, 1432-1451.
}
\keyword{datasets}
