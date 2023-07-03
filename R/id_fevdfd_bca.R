#' Identify...
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

  if (svar$type == "const") {
    A_hat <- svar$A_hat[, -1]
  } else {
    A_hat <- svar$A_hat
  }

  MY <- cbind(diag(k), matrix(0, k, k))
  MX <- rbind(A_hat, MY)
  ME <- rbind(svar$B, matrix(0, k, k))

  nx <- dim(MX)[[2]]

  gl <- 1000 # 1024
  freq_grid <- seq(0, 2 * pi, length.out = gl)
  freq_keep1 <- freq_grid >= min(freqs) & freq_grid <= max(freqs)
  freq_keep2 <- freq_grid >=  2 * pi - max(freqs) & freq_grid <= 2 * pi - min(freqs)
  freq_keep <- freq_keep1 | freq_keep2

  zi <- exp(-1i * freq_grid)
  r2pi <- 1 / (2 * pi)

  sp <- matrix(as.complex(0), gl, 1)
  sp2 <- matrix(as.complex(0), gl, k * k)

  for (gp in 1:gl) {
    if (freq_keep[gp] == 1) {
        fom <- t(MY[ti, ]) %*% (solve(diag(nx) - MX * zi[gp]) %*% ME)
        tmp <- r2pi * (fom %*% Conj(t(fom)))
        tmp <- freq_keep[gp] * tmp
        sp[gp] <- tmp
        tmp <- r2pi * (Conj(t(fom)) %*% fom)
        tmp <- freq_keep[gp] * tmp
        sp2[gp, ] <- Conj(t(c(tmp)))
    }
  }

    VTtmp <- 2 * pi * Re(stats::fft(sp, inverse = TRUE) / gl)
    VDtmp <- apply(sp2, 2, FUN = function(x) 2 * pi * Re(pracma::ifft(x)))

    contributions <- matrix(VDtmp[1, ] / VTtmp[1], k, k)

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
