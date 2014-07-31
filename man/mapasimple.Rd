\name{mapasimple}
\alias{mapasimple}
\title{MAPA Single Origin Forecast Wrapper}
\description{
The following function estimates MAPA and produces out-of-sample forecasts. 
It is not suggested, unless only simple single origin forecasts are required.
}
\usage{
mapasimple(y, ppy, fh = ppy, minimumAL = 1, maximumAL = ppy, comb = "mean",
	         output = "forecast", paral = 0, display = 0, outplot = 1, 
           hybrid = TRUE, model = "ZZZ")
}
\arguments{
  \item{y}{
    In sample observations of a time series (vector). If y == "paper" then 
	it prints paper reference.
}
  \item{ppy}{
    Periods in a season of the time series at the sampled frequency.
    If insample is a ts object then this is taken from its frequency, unless overriden.
}
  \item{fh}{
    Forecast horizon. Default = ppy.
}
  \item{minimumAL}{
    Lowest aggregation level to use. Default = 1, maximumAL>1.
}
  \item{maximumAL}{
    Highest aggregation level to use. Default = ppy.
}
  \item{comb}{
    Combination operator. One of "mean" or "median". Default is "mean"
}
  \item{output}{
    Type of output. One of "forecast" or "all". Default is "forecast". If output="all", 
	both forecasts and components estimates per aggregation level are provided.
}
  \item{paral}{
    Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 2 = yes 
	and initialise cluster. Default is 0.
}
  \item{display}{
    Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
}
  \item{outplot}{
    Provide output plot. 0 = no; 1 = time series and forecast only; 2 = time series, forecasts and components. For the components the rainbow colouring scheme is used. Red is aggregation level 1, followed by yellow, green, cyan, blue and magenta for the higher aggregation levels. Default is 1. 
}
  \item{hybrid}{
    Provide hybrid forecasts, as in Kourentzes et al. paper. If minimumAL > 1 then the minimumAL ETS forecasts are used. Default is TRUE.
}
  \item{model}{
    Allow only that type of ETS at each aggregation level. This follows similar coding to the ets function. The first letter refers to the error type ("A", "M" or "Z"); the second letter refers to the trend type ("N","A","Ad","M","Md" or "Z"); and the third letter refers to the season type ("N","A","M" or "Z"). The letters mean: "N"=none, "A"=additive, "M"=multiplicative and "Z"=automatically selected. A "d" for trend implies damped. By default model="ZZZ". If due to sample limitation ETS cannot be calculated at an aggregation level for the selected model, then no estimation is done for that specific level. For aggregation levels that seasonality becomes 1 then a non-seasonal model is estimated.
}
}
\value{
\item{forecasts}{Vector with forecasts.}
\item{components}{Array with MAPA components, if output="all".}
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
\code{\link{mapa}}.
}
\examples{
mapasimple(admissions)
}
\keyword{ ~mapa }
