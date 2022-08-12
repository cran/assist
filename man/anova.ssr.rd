\name{anova.ssr}
\alias{anova.ssr}
\title{
Testing a Non-parametric Function Fitted via Smoothing Splines
}
\description{
For smoothing spline models with a single smoothing parameter, test
the hypothesis that the unknown funciton 
lies in the null space using the local most powerful (LMP) test and 
a GCV or GML test.
}
\usage{
\method{anova}{ssr}(object, simu.size=100, ...)
}
\arguments{
\item{object}{
an object of class "ssr" fitted with a single smooting parameter.
}
\item{ simu.size}{
   an optional integer giving the number of simulations to calcualte 
p-values based on simulation. Default is 100.
  }
\item{\dots}{other available arguments, currently unused.}
}

\details{
For Gaussian data with one smoothing parameter, test the hypothesis that the function is in the 
null space \eqn{H_0}, i.e. the parametric part of the fitted model is sufficient. 
Available are the LMP and GCV or GML methods. However, the p-values cannot be calculated analytically 
since the null distributions for these testing statistics under \eqn{H_0} are unknown. 
Monte Carlo simulation is used to approximate the p-values
for the LMP, and GCV (if spar="v") or GML (if spar="m") methods. Due to computation burden, 
the smoothing parameters are fixed at their estimate in the currect calculation. 

When spar="m", an approximate p-value based on a mixture of two Chi-square distributions is also provided for the GML test, 
which tends to be conservationve (Pinherio and Bates, 2002). 

Methods further developed in  Liu and Wang (2004) and Liu, Meiring and Wang (2004) will be implemented in the future.

}
\value{
   a list including test values.
}
\author{ Chunlei Ke \email{chunlei_ke@yahoo.com} and Yuedong Wang \email{yuedong@pstat.ucsb.edu}}
\references{
Cox, D. and Koh, E. (1989). A smoothing spline based test of model adequency
in polynomial regression. Ann. Ins. Stat. Math. 41, 383-400.

Cox, D., Koh, E., Wahba, G. and Yandell, B.S. (1988). Testing the parameteric null
model hypothesis in semi-parametric partial and generalized spline models. 
Ann. Statist. 16, 113-119.

Wahba, G. (1990). Spline Models for Observational Data. SIAM, Vol. 59.

Pinherio, J. C. and Bates, D. M. (2000) Mixed-effects Models in S and S-Plus. Springer.

Liu, A. and Wang, Y. (2004) Hypothesis Testing in Smoothing Spline Models. 
Journal of Statistical Computation and Simulation, to appear. 

Liu, A., Meiring, W. and Wang, Y. (2004), Testing Generalized Linear Models Using Smoothing Spline Methods. Statistica Sinica, to appear,
}
\seealso{
\code{\link{ssr}}, \code{\link{print.anova.ssr}}
}
\examples{
\dontrun{
data(acid)

# fit a partial thin-plate spline
temp <- ssr(ph~t1+x1+x2, rk=tp(t1), data=acid, spar="m")
anova(temp, 500)
}
}
\keyword{file}
