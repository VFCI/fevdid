
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fevdid

<!-- badges: start -->

[![R-CMD-check](https://github.com/vfci/fevdid/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vfci/fevdid/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/VFCI/fevdid/branch/main/graph/badge.svg?token=2UCPW70685)](https://codecov.io/gh/VFCI/fevdid)
<!-- badges: end -->

R package to identify structural VAR shocks using maximization of
explained forecast error variances. Implemented to target either the
time domain or frequency domain.

## Installation

You can install the development version of fevdid from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("VFCI/fevdid")
```

## Usage

This is a basic example which shows you how to solve a common problem:

``` r
library(fevdid)
## basic example code
```

## References

[Chapter 6: “Forecast error variance
maximization”](https://jrenne.github.io/IdentifStructShocks/forecast-error-variance-maximization.html)
of the book *The Identification of Dynamic Stuctural Shocks*, by
Jean-Paul Renne and Kenza Benhima.

- [jrenne Github](https://github.com/jrenne)

Explains and codes up an example of identifying a shock that maximizes
forecast error variance contribution.
