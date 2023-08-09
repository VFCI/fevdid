#' Translation from the Business Cycle Anatomy replication code,
#' used to verify the output of id_fevdtd.
#'
#' @noRd
#' 
#' @param var, vars::VAR object
#' @param target, variable name or index to maximize its fevd
#' @param horizon, integer vector (can be length 1) of the horizon to maximize
#'
#' @return structural var
#'
id_fevdtd_bca <- function(var, target, horizon) {
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

  if (!is.numeric(horizon)) stop("Please provide an integer valued horizon.")
  if (!all(horizon > 0)) stop("Please provide only positive horizon values.")

  ## Fit a Choleskey SVAR (need orthogonal shocks)
  svar <- svars::id.chol(var)
  svar$B <- t(chol(stats::cov(stats::residuals(var))))

  ## Calculate IRFs out to horizon (then adj to 3-dim matrix from DF)
  irf <- vars::irf(svar, n.ahead = max(horizon))[[1]][, -1] |>
    apply(1, matrix, simplify = FALSE, nrow = k, ncol = k, byrow = TRUE) |>
    simplify2array()

  vtmp <- t(irf[ti, , ])

  for (j in min(horizon):max(horizon)) {
    V0 <- matrix(0, k, k)
    V1 <- 0
    for (i in 1:j) {
      V0 <- V0 + vtmp[i, ] %*% t(vtmp[i, ])
      V1 <- V1 + t(vtmp[i, ]) %*% vtmp[i, ]
    }
  }

  contributions <- V0 / V1[[1]]

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
  mvar$method <- "fevdtd"

  return(mvar)
}

#' Translation from the Business Cycle Anatomy replication code,
#' used to verify the output of id_fevdfd.
#'
#' @noRd
#' 
#' @param var, vars::VAR object
#' @param target, variable name or index to maximize its fevd
#' @param freqs vector of length 2 of min and max frequencies (0:2pi)
#'
#' @return structural var
#'
id_fevdfd_bca <- function(var, target, freqs) {
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
  svar1 <- as_statespace_var(svar$A_hat, svar$B)

  nx <- dim(svar1$mx)[[2]]

  gl <- 1000 # 1024 in original code

  ## Create grid of frequencies, set to True those to target
  freq_grid <- seq(0, 2 * pi, length.out = gl)
  f1 <- freq_grid >= min(freqs) & freq_grid <= max(freqs)
  f2 <- freq_grid >= 2 * pi - max(freqs) & freq_grid <= 2 * pi - min(freqs)
  freq_keep <- f1 | f2

  zi <- exp(-1i * freq_grid)
  r2pi <- 1 / (2 * pi)

  sp <- matrix(as.complex(0), gl, 1)
  sp2 <- matrix(as.complex(0), gl, k * k)

  for (gp in 1:gl) {
    if (freq_keep[gp] == 1) {
      fom <-
        t(svar1$my[ti, ]) %*% (solve(diag(nx) - svar1$mx * zi[gp]) %*% svar1$me)
      tmp <- r2pi * (fom %*% Conj(t(fom)))
      tmp <- freq_keep[gp] * tmp
      sp[gp] <- tmp
      tmp <- r2pi * (Conj(t(fom)) %*% fom)
      tmp <- freq_keep[gp] * tmp
      sp2[gp, ] <- Conj(t(c(tmp)))
    }
  }

  vt <- 2 * pi * Re(stats::fft(sp, inverse = TRUE) / gl)
  vd <- apply(sp2, 2, FUN = function(x) 2 * pi * Re(pracma::ifft(x)))

  contributions <- matrix(vd[1, ] / vt[1], k, k)

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

#' Identify the main shock by targetting the forecast error
#' variance contribution in the frequency domain.
#' Approximates the frequency fev by iterating out a length of
#' 'hmax' in time domain first.
#'
#' @noRd
#' 
#' @param var, vars::VAR object
#' @param target, variable name or index to maximize its fevd
#' @param freqs vector of length 2 of min and max frequencies (0:2pi)
#' @param hmax length of irfs to calculate, longer for better approximation
#'
#' @return main shock var
#'
id_fevdfd_approx <- function(var, target, freqs, hmax = 1000) {
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

  ## Calculate IRFs out to horizon (then adj to 3-dim matrix from DF)
  irf <- vars::irf(svar, n.ahead = hmax)[[1]][, -1] |>
    apply(1, matrix, simplify = FALSE, nrow = k, ncol = k, byrow = TRUE) |>
    simplify2array()

  ## Target matrix
  tm <- matrix(0, k, k)
  tm[ti, ti] <- 1

  ## Squared IRF contributions
  irf2 <- array(0, dim = c(k, k, hmax))
  irf2[, , 1] <- t(irf[, , 1]) %*% tm %*% irf[, , 1]
  for (h in 2:hmax) {
    irf2[, , h] <- t(irf[, , h]) %*% tm %*% irf[, , h]
  }

  freq_grid <- seq(0, 2 * pi, length.out = hmax)
  f1 <- freq_grid >= min(freqs) & freq_grid <= max(freqs)
  f2 <- freq_grid >= 2 * pi - max(freqs) & freq_grid <= 2 * pi - min(freqs)
  freq_keep <- f1 | f2

  ## Convert to Freq Domain
  fd_irf <- array(0, dim = c(k, k, hmax))
  for (i in 1:k) {
    for (j in 1:k) {
      td_vals <- irf[i, j, ]

      fd_vals <- stats::fft(td_vals)

      fd_keep <- fd_vals * as.integer(freq_keep)

      fd_irf[i, j, ] <- fd_keep
    }
  }

  ## Get squareds
  fd_irf_sq <- array(0, dim = c(k, k, hmax))
  for (t in 1:hmax) {
    fd_irf_sq[, , t] <- fd_irf[ti, , t] %*% Conj(t(fd_irf[ti, , t]))
  }

  contributions <- matrix(0, k, k)
  for (t in 1:hmax) {
    contributions <- contributions + Re(fd_irf_sq[, , t])
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
  mvar$method <- "id_fevdfd"
  mvar$target <- target
  mvar$freqs <- freqs

  class(mvar) <- c("fevdvar", "svars")

  return(mvar)
}
