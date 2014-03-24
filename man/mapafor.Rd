\name{mapafor}
\alias{mapafor}
\title{MAPA Forecast}
\description{
The following function produces in- and out-of-sample MAPA forecasts, for multiple steps ahead.
 This is the recommended function to use in forecasting with MAPA.
}
\usage{
mapafor(insample, mapafit, fh = -1, ifh = 1, comb = "mean", outplot = 1, hybrid = TRUE)
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
  \item{ifh}{
    In-sample forecast horizon. Default = 0.
}
  \item{comb}{
    Combination operator. One of "mean" or "median". Default is "mean".
}
  \item{outplot}{
    Provide output plot. 0 = no; 1 = yes. Default is 1. 
}
  \item{hybrid}{
    Provide hybrid forecasts, as in Kourentzes et al. paper. Default is TRUE.
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
\examples{
mapafit <- mapaest(admissions,outplot=0)
out <- mapafor(admissions,mapafit)
}
\keyword{ ~mapa }
\keyword{ ~mapacalc }
\keyword{ ~mapaest }
