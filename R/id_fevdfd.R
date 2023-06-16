#' Identify...
#'
#' @param var, vars::VAR object
#' @param target, variable name or index to maximize its fevd
#' @param horizon, integer vector (can be length 1) of the horizon to maximize
#'
#' @return structural var
#' @export
#'
#' @examples
#' x <- svars::USA
#' v <- vars::VAR(x, p = 2)
#' mvar <- id_fevdtd(v, "pi", 4:10)
#'
id_fevdfd <- function(var, target, freqs, hmax = 1000) {
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
  svar$B <- t(chol(cov(residuals(var))))

  ## Calculate IRFs out to horizon (then adj to 3-dim matrix from DF)
  irf <- vars::irf(svar, n.ahead = hmax)[[1]][, -1] |>
    apply(1, matrix, simplify = FALSE, nrow = k, ncol = k, byrow = TRUE) |>
    simplify2array()

  ## Target matrix
  tm <- matrix(0, k, k)
  tm[ti, ti] <- 1

  ## Squared IRF contributions
  irf2 <- array(0, dim = c(k, k, hmax))
  for (h in 1:hmax) {
    irf2[, , h] <- t(irf[, , h]) %*% tm %*% irf[, , h]
  }

  freq_grid <- 2 * pi * (0:(hmax - 1)) / hmax
  freq_keep1 <- freq_grid >= min(freqs) & freq_grid <= max(freqs)
  freq_keep2 <- freq_grid >=  2 * pi - max(freqs) & freq_grid <= 2 * pi - min(freqs)
  freq_keep <- freq_keep1 | freq_keep2

  contributions <- matrix(0, k, k)
  for (i in 1:k) {
    for (j in 1:k) {
        td_vals <- irf2[i, j, ]

        fd_vals <- fft(td_vals)

        fd_keep <- fd_vals * as.integer(freq_keep)

        td_keep <- Re(pracma::ifft(fd_keep))

        contributions[i, j] <- td_keep[1]
    }
  }

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
