#' Calculate the forecast error variance in the frequency domain.
#'
#' @param var, vars::VAR object
#' @param freqs vector of length 2 of min and max frequencies (0:2pi)
#' @param grid_size how fine the grid to approximate the frequency domain
#' @param fev Boolean true to return fev not fevd
#'
#' @return Matrix of forecast error variance decomposition in frequency domain
#' Indexed: frequencies, variables, shocks
#' @export
#'
#' @examples
#' x <- svars::USA
#' v <- vars::VAR(x, p = 2)
#' sv <- svars::id.chol(v)
#' fevd <- fevdfd(sv)
#'
fevdfd <- function(var, freqs = c(0, 2 * pi), grid_size = 1000, fev = FALSE) {
  ## Check parameter values are what is expected
  if (!inherits(var, "svars")) stop("Please pass a SVAR.")

  n <- colnames(var$y)
  k <- length(n)

  if (!is.numeric(freqs)) stop("Please provide numeric freqs.")
  if (!all(freqs >= 0 & freqs <= 2 * pi)) {
    stop("Please provide freqs between 0 and 2pi.")
  }

  ## Contstruct VAR(1) objects
  svar1 <- svar_to_svar1(var)

  nx <- dim(svar1$mx)[[2]]

  gl <- grid_size

  ## Create grid of frequencies, set to True those to target
  freq_grid <- seq(0, 2 * pi, length.out = gl)
  f1 <- freq_grid >= min(freqs) & freq_grid <= max(freqs)
  f2 <- freq_grid >= 2 * pi - max(freqs) & freq_grid <= 2 * pi - min(freqs)
  freq_keep <- f1 | f2

  zi <- exp(-1i * freq_grid)
  r2pi <- 1 / (2 * pi)

  freq_fev <- array(0, dim = c(gl, k, k))
  freq_totvar <- array(0, dim = c(gl, k))

  for (gp in 1:gl) {
    if (freq_keep[gp] == 1) {
      x2 <-
        svar1$my %*% (solve(diag(nx) - svar1$mx * zi[gp]) %*% svar1$me)
      x_sq2 <- r2pi * abs(x2)^2
      freq_fev[gp, , ] <- freq_keep[gp] * x_sq2

      freq_totvar[gp, ] <- rowSums(Re(x_sq2))
    }
  }

  if (fev == TRUE) {
    return(freq_fev)
  }

  freq_fevd <- array(0, dim = c(gl, k, k))
  for (gp in 1:gl) {
    freq_fevd[gp, , ] <- freq_fev[gp, , ] / freq_totvar[gp, ]
  }

  return(freq_fevd)
}
