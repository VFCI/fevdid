## Install bcadata
# devtools::install_github("VFCI/bcadata")
require(dplyr)

## Load in the BCA data
## Turns out they only use up to 2017:Q1, for the classical VAR,
## for the Bayesian ones they use through 2017:Q4
x <- bcadata::original_bcadata |>
  filter(date <= as.Date("2017-01-01")) |>
  select(-date)

## Fit VAR(2), find Choleskey for comparison
v <- vars::VAR(x, 2, type = "const")
sv <- svars::id.chol(v)


## Main Shock, target unemployment, Time Domain
## Their code, but translated to R
mtv1 <- id_fevdtd_bca(v, "unemployment", 1:5)

## My code, gives different values
mtv2 <- id_fevdtd(v, "unemployment", 1:5)

## My code, gives same value
mtv3 <- id_fevdtd(v, "unemployment", 5)

## The original BCA code seems to have an error where it incorrectly
## targets the max horizon only, not the range.
print(mtv1); print(mtv2); print(mtv3)


## Main shock, target unemployment, Freq Domain
bc_freq <- c(2 * pi / 32, 2 * pi / 6)
mfv1 <- id_fevdfd(v, "unemployment", bc_freq)
mfv2 <- id_fevdfd_bca(v, "unemployment", bc_freq)

## Different, but similar.
## Need to understand what their code is doing better.
print(mfv1); print(mfv2)
