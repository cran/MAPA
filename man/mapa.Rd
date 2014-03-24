\name{mapa}
\alias{mapa}

\title{Multiple Aggregation Prediction Algorithm (Wrapper)}

\description{Wrapper to estimate and produce MAPA in- and out-of-sample forecasts. 
Uses mapaest and mapafor.}

\usage{
mapa(insample, ppy, fh = ppy, ifh = 1, minimumAL = 1, maximumAL = ppy, 
	comb = "mean", paral = 0, display = 0, outplot = 1, hybrid = TRUE, model = "ZZZ")
}

\arguments{
  \item{insample}{
    In sample observations of a time series (vector). If insample == "paper" then it prints paper reference.
}
  \item{ppy}{
    Periods in a season of the time series at the sampled frequency.
    If insample is a ts object then this is taken from its frequency, unless overriden.
}
  \item{fh}{
    Forecast horizon. Default = ppy.
}
  \item{ifh}{
    Lower aggregation level to use. Default = 1.
}
  \item{minimumAL}{
    Lower aggregation level to use. Default = 1.
}
  \item{maximumAL}{
    Highest aggregation level to use. Default = ppy, maximumAL>1.
}
  \item{comb}{
    Combination operator. One of "mean" or "median". Default is "mean".
}
  \item{paral}{
    Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 2 = yes and initialise cluster. Default is 0.
}
  \item{display}{
    Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
}
  \item{outplot}{
    Provide output plot. 0 = no; 1 = yes. Default is 1.
}
  \item{hybrid}{
    Provide hybrid forecasts, as in Kourentzes et al. paper. If minimumAL > 1 then the minimumAL ETS forecasts are used. Default is TRUE.
}
  \item{model}{
    Allow only that type of ETS at each aggregation level. This follows similar coding to the ets function. The first letter refers to the error type ("A", "M" or "Z"); the second letter refers to the trend type ("N","A","Ad","M","Md" or "Z"); and the third letter refers to the season type ("N","A","M" or "Z"). The letters mean: "N"=none, "A"=additive, "M"=multiplicative and "Z"=automatically selected. A "d" for trend implies damped. By default model="ZZZ". If due to sample limitation ETS cannot be calculated at an aggregation level for the selected model, then no estimation is done for that specific level. 
}
}
\value{
\item{out$infor }{In-sample forecasts}
\item{out$outfor }{Out-of-sample forecasts}
\item{out$MSE }{In-sample MSE error}
\item{out$MAE }{In-sample MAE error}
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
\code{\link{mapaest}}, \code{\link{mapafor}}, \code{\link{mapasimple}}.
}
\examples{
out <- mapa(admissions)
}
\keyword{ ~mapaest }
\keyword{ ~mapafor }
\keyword{ ~mapasimple }