#' Calculate either the IRF in the frequency domain.
#'
#' @param var A VAR. Currently supports either a 'svars' or 'fevdvar' object.
#' @param ... Not currently used.
#'
#' @return irffd object, storing a data.frame
#' @export
#'
#' @examples
#' x <- svars::USA
#' v <- vars::VAR(x, p = 2)
#' sv <- svars::id.chol(v)
#' irffd <- irffd(sv)
#'
irffd <- function(var, ...) {
  UseMethod("irffd")
}


#' Find the IRF in FD for a VAR
#'
#' @param var A VAR. Currently supports either a 'svars' or 'fevdvar' object.
#' @param freqs vector of length 2 of min and max frequencies (0:2pi)
#' @param grid_size how fine the grid to approximate the frequency domain
#'
#' @return list with vector of frequencies and matrix of irf values
#'
irffd_generic <- function(
    var,
    freqs = c(0, 2 * pi),
    grid_size = 1000) {
  n <- colnames(var$y)
  k <- length(n)

  impulse_names <- n
  response_names <- n

  ## Set as factors
  impulse_names <-
    factor(impulse_names, levels = impulse_names, ordered = TRUE)
  response_names <-
    factor(response_names, levels = response_names, ordered = TRUE)

  if (!is.numeric(freqs)) stop("Please provide numeric freqs.")
  if (!all(freqs >= 0 & freqs <= 2 * pi)) {
    stop("Please provide freqs between 0 and 2pi.")
  }

  ## Contstruct VAR(1) objects
  svar1 <- as_statespace_var(var$A_hat, var$B)

  nx <- dim(svar1$mx)[[2]]

  gl <- grid_size

  ## Create grid of frequencies, set to True those to target
  freq_grid <- seq(0, 2 * pi, length.out = gl)
  f1 <- freq_grid >= min(freqs) & freq_grid <= max(freqs)
  f2 <- freq_grid >= 2 * pi - max(freqs) & freq_grid <= 2 * pi - min(freqs)
  freq_keep <- f1 | f2

  zi <- exp(-1i * freq_grid)
  r2pi <- 1 / (2 * pi)

  freq_irf <- array(0, dim = c(gl, k, k))

  for (gp in 1:gl) {
    if (freq_keep[gp] == 1) {
      irf <-
        svar1$my %*% (solve(diag(nx) - svar1$mx * zi[gp]) %*% svar1$me)
      irf <- r2pi * Re(irf)
      freq_irf[gp, , ] <- freq_keep[gp] * irf
    }
  }

  return(list(
    freq_grid = freq_grid,
    freq_irf = freq_irf
  ))
}


#' Method to calculate irffd for fevdvar (this package)
#'
#' @param freqs vector of length 2 of min and max frequencies (0:2pi)
#' @param grid_size how fine the grid to approximate the frequency domain
#'
#' @rdname irffd
#' @name irffd
#' @aliases irffd.fevdvar
#' @export
#'
irffd.fevdvar <- function(
    var,
    freqs = c(0, 2 * pi),
    grid_size = 1000,
    ...) {
  k <- var$K

  irffd <- irffd_generic(var, freqs, grid_size)

  response_names <- colnames(var$y)
  impulse_names <- var$impulse_names

  df <- data.frame(
    f = rep(irffd$freq_grid, times = k * k),
    impulse = rep(impulse_names, each = grid_size * k),
    response = rep(response_names, each = grid_size, times = k),
    irf = c(irffd$freq_irf)
  )

  irffd <- list(irffd = df)
  class(irffd) <- "irffd"

  return(irffd)
}


#' Method to calculate irffd for svars (id.chol)
#'
#' @param freqs vector of length 2 of min and max frequencies (0:2pi)
#' @param grid_size how fine the grid to approximate the frequency domain
#'
#' @rdname irffd
#' @name irffd
#' @aliases irffd.svars
#' @export
#'
irffd.svars <- function(
    var,
    freqs = c(0, 2 * pi),
    grid_size = 1000,
    ...) {
  k <- var$K

  irffd <- irffd_generic(var, freqs, grid_size)

  response_names <- colnames(var$y)
  impulse_names <- colnames(var$y)

  df <- data.frame(
    f = rep(irffd$freq_grid, times = k * k),
    impulse = rep(impulse_names, each = grid_size * k),
    response = rep(response_names, each = grid_size, times = k),
    irf = c(irffd$freq_irf)
  )

  irffd <- list(irffd = df)
  class(irffd) <- "irffd"

  return(irffd)
}
