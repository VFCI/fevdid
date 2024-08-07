% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forecast.R
\name{forecast}
\alias{forecast}
\title{Calculate the forecasts for a VAR out to a specified horizon
Conditions on information at time "t"}
\usage{
forecast(var, horizon = 1)
}
\arguments{
\item{var}{vars::VAR or svars object}

\item{horizon}{number of steps out to calculate forecast}
}
\value{
forecast
}
\description{
Calculate the forecasts for a VAR out to a specified horizon
Conditions on information at time "t"
}
\examples{
x <- svars::USA
v <- vars::VAR(x, p = 2)
forecast(v, 10)

}
