#' Identify the main shock by targetting the forecast error
#' variance contribution in the frequency domain.
#' 
#' @param var, vars::VAR object
#' @param target, variable name or index to maximize its fev
#' @param freqs vector of length 2 of min and max frequencies (0:2pi)
#' @param grid_size how fine the grid to approximate the frequency domain
#'
#' @return structural var
#' @export
#'
#' @examples
#' x <- svars::USA
#' v <- vars::VAR(x, p = 2)
#' mvar <- id_fevdfd(v, "pi", c(2 * pi / 32, 2 * pi / 6))
#'
id_fevdfd <- function(var, target, freqs, grid_size = 1000) {
  ## Check parameter values are what is expected
  if (!inherits(var, "varest")) stop("Please pass a VAR from 'vars::VAR'.")

  n <- colnames(var$y)
  k <- length(n)
  ni <- 1:k

  if (is.numeric(target)) {
    if (!target %in% ni) stop("Please provide a valid target variable.")
    ti <- target
    t <- n[target]
  } else {
    target <- as.character(target)
    if (!target %in% n) stop("Please provide a valid target variable.")
    ti <- which(n %in% target)
    t <- target
  }

  if (!is.numeric(freqs)) stop("Please provide numeric freqs.")
  if (!all(freqs > 0 & freqs < 2 * pi)) {
    stop("Please provide freqs between 0 and 2pi.")
  }

  ## Fit a Choleskey SVAR (need orthogonal shocks)
  svar <- svars::id.chol(var)
  svar$B <- t(chol(stats::cov(stats::residuals(var))))

  ## Contstruct VAR(1) objects
  svar1 <- svar_to_svar1(svar)

  nx <- dim(svar1$mx)[[2]]

  gl <- grid_size

  ## Create grid of frequencies, set to True those to target
  freq_grid <- seq(0, 2 * pi, length.out = gl)
  f1 <- freq_grid >= min(freqs) & freq_grid <= max(freqs)
  f2 <- freq_grid >= 2 * pi - max(freqs) & freq_grid <= 2 * pi - min(freqs)
  freq_keep <- f1 | f2

  zi <- exp(-1i * freq_grid)
  r2pi <- 1 / (2 * pi)

  freq_fev <- matrix(as.complex(0), gl, k * k)

  for (gp in 1:gl) {
    if (freq_keep[gp] == 1) {
      x <-
        t(svar1$my[ti, ]) %*% (solve(diag(nx) - svar1$mx * zi[gp]) %*% svar1$me)
      x_sq <- r2pi * Conj(t(x)) %*% x
      x_sq <- freq_keep[gp] * x_sq
      freq_fev[gp, ] <- Conj(t(c(x_sq)))
    }
  }

  contributions <- matrix(colSums(Re(freq_fev)), k, k)

  ## Max eigen value
  e <- eigen(contributions)
  ei <- which(e$value == max(e$value))
  evec <- e$vectors[, ei]

  ## Construct max rotation matrix
  q <- matrix(0, k, k)
  q[, 1] <- evec
  q[, 2:k] <- pracma::nullspace(t(evec))

  ## Insert resulting matrix into var
  mvar <- svar
  mvar$B <- svar$B %*% q
  mvar$method <- "fevdfd"

  return(mvar)
}
