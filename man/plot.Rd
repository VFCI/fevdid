% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot}
\alias{plot}
\alias{plot.fevdfd}
\alias{plot.fevdvarfevd}
\title{Plot the forecast variance decomposition.}
\usage{
\method{plot}{fevdfd}(x, y, stacked = TRUE, vlines = NULL, ...)

\method{plot}{fevdvarfevd}(x, y, stacked = TRUE, vlines = NULL, ...)
}
\arguments{
\item{x}{object of class "svarfevd"}

\item{y}{Not used.}

\item{stacked}{Boolean.
True for stacked columns, False (default) for unstacked line chart.}

\item{vlines}{Vector of x-axis points at which to draw a vline.
Useful for highlighting areas.}

\item{...}{Currently not used.}
}
\value{
ggplot of fevdfd

ggplot of fevd
}
\description{
Alias for the function plot in frequency domain.
}
