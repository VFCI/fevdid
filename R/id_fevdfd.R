#' Identify the main shock by targetting the forecast error
#' variance contribution in the frequency domain.
#'
#' @param x either a `vars::VARS` or a `BVAR::bvar` object
#' @param target, variable name or index to maximize its fev
#' @param freqs vector of length 2 of min and max frequencies (0:2pi)
#' @param grid_size how fine the grid to approximate the frequency domain
#' @param freq_grid Default to NULL.
#' Pass a vector of frequencies in c(0,2pi) to target those specific
#' frequencies.  Overides `freqs` and `grid_size` arguments.
#' @param sign Default to "positive". Can be "negative".  Ensures the
#' cummulative impact of the main shock on the target variable is the
#' given sign.
#' @param sign_horizon Default to 1. The horizon through which to accumulate the
#' impact of the shock.
#'
#' @return structural var
#' @export
#'
#' @examples
#' x <- svars::USA
#' v <- vars::VAR(x, p = 2)
#' mvar <- id_fevdfd(v, "pi", c(2 * pi / 32, 2 * pi / 6))
#'
id_fevdfd <- function(
  x,
  target,
  freqs = c(0, 2 * pi),
  grid_size = 1000,
  freq_grid = NULL,
  sign = "positive",
  sign_horizon = 1
) {
  UseMethod("id_fevdfd")
}


#' @rdname id_fevdfd
#' @name id_fevdfd
#' @aliases id_fevdfd.varest
#'
#' @export
id_fevdfd.varest <- function(
  x,
  target,
  freqs = c(0, 2 * pi),
  grid_size = 1000,
  freq_grid = NULL,
  sign = "positive",
  sign_horizon = 1
) {
  n <- colnames(x$y)
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

  if (!is.null(freq_grid)) {
    if (!is.numeric(freqs)) stop("Please provide numeric freqs.")
    if (!all(freqs >= 0 & freqs <= 2 * pi)) {
      stop("Please provide freqs between 0 and 2pi.")
    }
  }

  ## Fit a Choleskey SVAR (need orthogonal shocks)
  svar <- svars::id.chol(x)
  sigma <- stats::cov(stats::residuals(x))
  betas <- svar$A_hat
  svar$B <- t(chol(sigma))

  q <- id_fevdfd_findq(betas, svar$B, ti, freqs, grid_size, freq_grid)

  ## Insert resulting matrix into var
  mvar <- svar
  mvar$Q <- q
  mvar$B <- svar$B %*% q
  mvar$method <- "id_fevdfd"
  mvar$target <- target
  mvar$freqs <- freqs

  class(mvar) <- c("fevdvar", "svars")

  ## Insure the sign is as expected
  irf <-
    vars::irf(mvar, impulse = "Main", response = t, n.ahead = sign_horizon)
  if (sign == "positive" || sign == "pos") {
    if (sum(irf$irf[1:sign_horizon, "irf"]) < 0) {
      mvar$Q <- -1 * mvar$Q
      mvar$B <- -1 * mvar$B
    }
  } else if (sign == "negative" || sign == "neg") {
    if (sum(irf$irf[1:sign_horizon, "irf"]) > 0) {
      mvar$Q <- -1 * mvar$Q
      mvar$B <- -1 * mvar$B
    }
  }

  return(mvar)
}


#' @rdname id_fevdfd
#' @name id_fevdfd
#' @aliases id_fevdfd.varboot
#'
#' @export
id_fevdfd.varboot <- function(
  x,
  target,
  freqs = c(0, 2 * pi),
  grid_size = 1000,
  freq_grid = NULL,
  sign = "positive",
  sign_horizon = 1
) {
  id_fevdfd.varest(x, target, freqs, grid_size, freq_grid, sign, sign_horizon)
}


#' @rdname id_fevdfd
#' @name id_fevdfd
#' @aliases id_fevdfd.bvartools
#'
#' @export
id_fevdfd.bvartools <- function(
  x,
  target,
  freqs = c(0, 2 * pi),
  grid_size = 1000,
  freq_grid = NULL,
  sign = "positive",
  sign_horizon = 1
) {
  k <- nrow(x$y)
  p <- ncol(x$A[, 1]) %% k
  var_names <- rownames(x$y)
  iterations <- ncol(x$Sigma)

  ti <- which(var_names == target)

  rots_sigma <- matrix(NA, k^2, iterations)

  ## Find the main shock rotation matrices, q, and rotate the sigma matrix
  for (i in 1:iterations) {
    if (i %% 10^2 == 0) print(i)
    a <- matrix(x$A[, i], k)
    c <- matrix(x$C[, i], k, 1)

    a_hat <- cbind(c, a)

    sigma <- matrix(x$Sigma[, i], k, k)
    b <- t(chol(sigma))

    q <- id_fevdfd_findq(a_hat, b, ti, freqs, grid_size, freq_grid)

    impulse <- b %*% q
    if (impulse[ti, 1] < 0) impulse <- impulse * -1

    rots_sigma[, i] <- c(impulse)
  }

  mvar <- x
  mvar$method <- "id_fevdfd"
  mvar$Sigma <- rots_sigma

  return(mvar)
}



#' @rdname id_fevdfd
#' @name id_fevdfd
#' @aliases id_fevdfd.bvar
#'
#' @export
#'
id_fevdfd.bvar <- function(
  x,
  target,
  freqs = c(0, 2 * pi),
  grid_size = 1000,
  freq_grid = NULL,
  sign = "positive",
  sign_horizon = 1
) {
  n <- colnames(x$meta$Y)
  k <- length(n)
  ni <- 1:k
  draws <- x$meta$n_save

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

  if (!is.null(freq_grid)) {
    if (!is.numeric(freqs)) stop("Please provide numeric freqs.")
    if (!all(freqs >= 0 & freqs <= 2 * pi)) {
      stop("Please provide freqs between 0 and 2pi.")
    }
  }

  ## Copy BVAR object
  mvar <- x

  ## Progress Bar
  tick <- 0
  ticks <- 20
  pb <- paste0(c(rep("=", tick), rep("-", ticks - tick)), collapse = "")
  cat("[", pb, "]")

  ## Iterate over Bayesian draws
  for (i in 1:draws) {
    ## Progress bar updates
    if (i %% (draws / ticks) == 0) {
      tick <- tick + 1
      pb <-
        paste0(c(rep("=", tick), rep("-", ticks - tick)), collapse = "")
      cat("\r[", pb, "]")
    }

    betas <- t(x$beta[i, , ])
    sigma <- x$sigma[i, , ]
    b <- t(chol(sigma))

    q <- id_fevdfd_findq(betas, b, ti, freqs, grid_size, freq_grid)


    impact <- b %*% q
    if (sign == "positive" || sign == "pos") {
      if (impact[ti, 1] < 0) {
        q <- -1 * q
      }
    } else if (sign == "negative" || sign == "neg") {
      if (impact[ti, 1] > 0) {
        q <- -1 * q
      }
    }

    mvar$sigma[i, , ] <- b %*% q
  }
  ## Insert resulting matrix into var
  mvar$method <- "id_fevdfd"
  mvar$target <- target
  mvar$freqs <- freqs

  class(mvar) <- c("fevdbvar", "bvar")

  return(mvar)
}

#'
#' Find the rotation matrix Q that maximizes fevd over
#' the given frequencies for the target variable.
#'
#' @param betas matrix of coeffecients from VAR
#' @param sigma variance-covariance matrix from VAR
#' @param target_index index of variable to target
#' @param freqs vector of length 2 of min and max frequencies (0:pi)
#' @param grid_size how fine the grid to approximate the frequency domain
#' @param freq_grid Default to NULL.
#' Pass a vector of frequencies in c(0,2pi) to target those specific
#' frequencies.  Overides `freqs` and `grid_size` arguments.
#'
#' @return matrix q
#'
id_fevdfd_findq <- function(
  betas,
  sigma,
  target_index,
  freqs,
  grid_size,
  freq_grid = NULL
) {
  ## Construct Objects
  k <- nrow(sigma)
  ti <- target_index

  ## Contstruct VAR(1) objects
  ssv <- as_statespace_var(betas, sigma)

  nx <- dim(ssv$mx)[[2]]

  ## Create grid of frequencies, set to True those to target
  if (is.null(freq_grid)) {
    freq_grid <- seq(0, 2 * pi, length.out = grid_size)
    f1 <- freq_grid >= min(freqs) & freq_grid <= max(freqs)
    f2 <- freq_grid >= 2 * pi - max(freqs) & freq_grid <= 2 * pi - min(freqs)
    freq_keep <- f1 | f2
  } else {
    freq_grid <- freq_grid
    grid_size <- length(freq_grid)
    freq_keep <- rep(TRUE, grid_size)
  }

  zi <- exp(-1i * freq_grid)
  r2pi <- 1 / (2 * pi)

  freq_fev <- matrix(as.complex(0), grid_size, k * k)

  for (gp in 1:grid_size) {
    if (freq_keep[gp] == 1) {
      x <- t(ssv$my[ti, ]) %*% (solve(diag(nx) - ssv$mx * zi[gp]) %*% ssv$me)
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

  return(q)
}
