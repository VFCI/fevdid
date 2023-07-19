
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

### Main Shock in the Time Domain

Example of identifying the main shock along the time domain. Here the
target variable is inflation over four to ten quarters out. The figure
compares the identified main shock to the Choleskey shock for inflation,
which already explains most of the variation in this simple 3-variable
VAR. Clearly, however, the contribution from the Main shock is higher
over the targeted time period.

``` r
require(fevdid)

## US data on inflation (pi), output (x), and federal funds rate (i)
x <- svars::USA
v <- vars::VAR(x, p = 2)
sv <- svars::id.chol(v, order_k = c("x", "pi", "i"))

## Find shock that maximizes forecast error variance for inflation (pi)
## from 4 to 10 quarters out
mv <- id_fevdtd(v, target = "pi", horizon = 4:10)

## Comparing fevds
fevdsv <- vars::fevd(sv, n.ahead = 20)$pi
fevdmv <- vars::fevd(mv, n.ahead = 20)$pi


## Plotting
require(ggplot2)
cbind(h = 1:20, fevdsv, fevdmv) |>
  tidyr::pivot_longer(cols = !h) |>
  dplyr::filter(name %in% c("Main", "pi")) |>
  ggplot(aes(x = h, y = value, color = name)) +
  geom_vline(xintercept = c(4, 10)) +
  geom_line()
```

<img src="man/figures/README-example_td-1.png" width="100%" />

Getting impulse responses from the identified main shock VARs is handled
by the usual `irf` function.

``` r

irfsv <- vars::irf(sv, n.ahead = 20)
irfmv <- vars::irf(mv, n.ahead = 20)

cowplot::plot_grid(plot(irfsv), plot(irfmv), nrow = 1)
```

<img src="man/figures/README-exanmple_td_irf-1.png" width="100%" />

Plotting all of the forecast variance decompositions is equally
straightforward (above the code is used to compare two separate fevds).

``` r

fevdsv <- vars::fevd(sv, n.ahead = 20)
fevdmv <- vars::fevd(mv, n.ahead = 20)

cowplot::plot_grid(plot(fevdsv), plot(fevdmv), nrow = 1)
```

<img src="man/figures/README-exanmple_td_fevd-1.png" width="100%" />

### Main Shock in the Frequency Domain

Example of identifying the main shock along the frequency domain. Here
the target variable is inflation between $\frac{2\pi}{32}$ and
$\frac{2\pi}{6}$, the “business cycle” frequencies. The figure compares
the identified main shock to the Choleskey shock for inflation, which
already explains most of the variation in this simple 3-variable VAR.
The contribution from the Main shock is higher over the targeted
frequency period.

``` r
require(fevdid)

## US data on inflation (pi), output (x), and federal funds rate (i)
x <- svars::USA
v <- vars::VAR(x, p = 2)
sv <- svars::id.chol(v, order_k = c("x", "pi", "i"))

## Find shock that maximizes forecast error variance for inflation (pi)
## in the "business cycle" frequencies
bc_freqs <- c(2 * pi / 32, 2 * pi / 6)
mfv <- id_fevdfd(v, target = "pi", freqs = bc_freqs)

## Comparing fevds
fevdsv <- fevdfd(sv)$pi
fevdmv <- fevdfd(mfv)$pi

## Plotting
require(ggplot2)
dplyr::full_join(fevdsv, fevdmv, by = "f") |>
  dplyr::as_tibble() |>
  tidyr::pivot_longer(cols = !f) |>
  dplyr::filter(name %in% c("Main", "pi")) |>
  ggplot(aes(x = f, y = value, color = name)) +
  geom_vline(xintercept = bc_freqs) +
  geom_line()
```

<img src="man/figures/README-example_fd-1.png" width="100%" />

There is also a generic plot function for all of the frequency domain
fevds.

``` r

fevdsv <- fevdfd(sv)
fevdmv <- fevdfd(mv)

cowplot::plot_grid(plot(fevdsv), plot(fevdmv), nrow = 1)
```

<img src="man/figures/README-exanmple_fd_fevd-1.png" width="100%" />

## References

[Chapter 6: “Forecast error variance
maximization”](https://jrenne.github.io/IdentifStructShocks/forecast-error-variance-maximization.html)
of the book *The Identification of Dynamic Stuctural Shocks*, by
Jean-Paul Renne and Kenza Benhima.

- [jrenne Github](https://github.com/jrenne)

Explains and codes up an example of identifying a shock that maximizes
forecast error variance contribution.

“[What moves real
GNP?](http://fmwww.bc.edu/repec/esNAWM04/up.2923.1054309431.pdf)”.
Harlad Uhlig. (2003).
