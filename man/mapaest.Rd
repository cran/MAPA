\name{mapaest}
\alias{mapaest}
\title{MAPA Estimation}
\description{
The following function estimates MAPA and saves all fitted model.
}
\usage{
mapaest(insample, ppy, minimumAL = 1, maximumAL = ppy, paral = 0, 
	display = 0, outplot = 1)
}
\arguments{
  \item{insample}{
    In sample observations of a time series (vector).
}
  \item{ppy}{
    Periods in a season of the time series at the sampled frequency.
    If insample is a ts object then this is taken from its frequency, unless overriden.
}
  \item{minimumAL}{
    Lower aggregation level to use. Default = 1, maximumAL>1.
}
  \item{maximumAL}{
    Highest aggregation level to use. Default = ppy.
}
  \item{paral}{
    Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 2 = yes 
	and initialise cluster. Default is 0.
}
  \item{display}{
    Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
}
  \item{outplot}{
    Provide output plot. 0 = no; 1 = yes. Default is 1.
}
}
\value{
\item{mapafit }{Estimated MAPA model structure}
}
\references{
Kourentzes N., Petropoulos F., Trapero J.R., 2014. Improving forecasting by estimating time 
series structural components across multiple frequencies, International Journal of Forecasting,
 30(2), 291-302.
}
\author{
Nikolaos Kourentzes and Fotios Petropoulos
}
\examples{
mapafit <- mapaest(admissions)
}
\keyword{ ~mapa }
\keyword{ ~mapafor }
