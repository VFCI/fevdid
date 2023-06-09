% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fevd.R
\name{fevd}
\alias{fevd}
\alias{fevd.fevdvar}
\title{Forecast error variance decomposition.}
\usage{
\method{fevd}{fevdvar}(x, n.ahead = 10, ...)
}
\arguments{
\item{x}{SVAR object of class "fevdvar"}

\item{n.ahead}{Integer specifying the steps.}

\item{...}{Currently not used.}
}
\value{
Matrix of forecast error variance decomposition in frequency domain
Indexed: frequencies, variables, shocks
}
\description{
Alias for the function from svars::fevd.svars,
so that the shock names are correct.
}
\examples{
x <- svars::USA
v <- vars::VAR(x, p = 2)
mvar <- id_fevdtd(v, "pi", 4:10)
vars::fevd(mvar)

}
