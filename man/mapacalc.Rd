\name{mapacalc}
\alias{mapacalc}
\title{MAPA Internal Calculation of Forecasts}
\description{
The following function produces a simple MAPA forecast from a given origin. 
This is meant to be an internal function. Use mapafor instead.
}
\usage{
mapacalc(insample, mapafit, fh = 0, comb = "mean", output = "forecast", outplot = 0, 
	hybrid = TRUE)
}
\arguments{
  \item{insample}{
    In sample observations of a time series (vector).
}
  \item{mapafit}{
    Fitted MAPA model (from mapaest).
}
  \item{fh}{
    Forecast horizon. Default = ppy.
}
  \item{comb}{
    Combination operator. One of "mean" or "median". Default is "mean".
}
  \item{output}{
    Type of output. One of "forecast" or "all". Default is "forecast". If output="all", 
	both forecasts and components estimates per aggregation level are provided.
}
  \item{outplot}{
    Provide output plot. 0 = no; 1 = time series and forecast only; 2 = time series, 
	forecasts and components. Default is 1.
}
\item{hybrid}{
    Provide hybrid forecasts, as in Kourentzes et al. paper. If minimumAL > 1 then the minimumAL ETS forecasts are used. Default is TRUE.
}
}
\value{
\item{forecasts }{Vector with forecasts}
\item{components }{Array with MAPA components, if output="all"}
}
\references{
Kourentzes N., Petropoulos F., Trapero J.R., 2014. Improving forecasting by estimating time 
series structural components across multiple frequencies, International Journal of Forecasting,
 30(2), 291-302.
}
\author{
Nikolaos Kourentzes and Fotios Petropoulos
}
\seealso{
\code{\link{mapafor}}, \code{\link{mapa}}.
}
\examples{
mapafit <- mapaest(admissions,outplot=0)
mapacalc(admissions,mapafit,outplot=2)
}
\keyword{ ~mapafor }
