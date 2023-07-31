#' Identify the main shock by targetting the forecast error
#' variance contribution in the time domain.
#'
#' @param x vars::VAR object
#' @param target variable name or index to maximize its fevd
#' @param horizon integer vector (can be length 1) of the horizon to maximize
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
#' mvar <- id_fevdtd(v, "pi", 4:10)
#'
id_fevdtd <- function(
    x,
    target,
    horizon,
    sign = "positive",
    sign_horizon = 1) {
  UseMethod("id_fevdtd")
}

#' @rdname id_fevdtd
#' @name id_fevdtd
#' @aliases id_fevdtd.varest
#'
#' @export
id_fevdtd.varest <- function(
    x,
    target,
    horizon,
    sign = "positive",
    sign_horizon = 1) {
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

  if (!is.numeric(horizon)) stop("Please provide an integer valued horizon.")
  if (!all(horizon > 0)) stop("Please provide only positive horizon values.")

  ## Fit a Choleskey SVAR (need orthogonal shocks)
  b <- t(chol(stats::cov(stats::residuals(x))))
  svar <- svars::id.chol(x)

  q <- id_fevdtd_findq(svar$A_hat, b, ti, horizon)

  ## Insert resulting matrix into var
  mvar <- svar
  mvar$Q <- q
  mvar$B <- b %*% q
  mvar$method <- "id_fevdtd"
  mvar$target <- target
  mvar$horizon <- horizon

  class(mvar) <- c("fevdvar", "svars")

  ## Insure the sign is as expected
  mssv <- as_statespace_var(mvar$A_hat, mvar$B)
  irf <- irf_ssv(mssv, n_ahead = sign_horizon)[ti, 1, ]

  if (sign == "positive" || sign == "pos") {
    if (sum(irf) < 0) {
      mvar$Q <- -1 * mvar$Q
      mvar$B <- -1 * mvar$B
    }
  } else if (sign == "negative" || sign == "neg") {
    if (sum(irf) > 0) {
      mvar$Q <- -1 * mvar$Q
      mvar$B <- -1 * mvar$B
    }
  }

  return(mvar)
}

#' @rdname id_fevdtd
#' @name id_fevdtd
#' @aliases id_fevdtd.varboot
#'
#' @export
id_fevdtd.varboot <- function(
    x,
    target,
    horizon,
    sign = "positive",
    sign_horizon = 1) {
  id_fevdtd.varest(x, target, horizon, sign, sign_horizon)
}


#' @rdname id_fevdtd
#' @name id_fevdtd
#' @aliases id_fevdtd.bvar
#'
#' @export
#'
id_fevdtd.bvar <- function(
    x,
    target,
    horizon,
    sign = "positive",
    sign_horizon = 1) {
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

  if (!is.numeric(horizon)) stop("Please provide an integer valued horizon.")
  if (!all(horizon > 0)) stop("Please provide only positive horizon values.")

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

    q <- id_fevdtd_findq(betas, b, ti, horizon)

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
  mvar$method <- "id_fevdtd"
  mvar$target <- target
  mvar$horizon <- horizon

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
#' @param horizon integer vector (can be length 1) of the horizon to maximize
#'
#' @return matrix q
#'
id_fevdtd_findq <- function(
    betas,
    sigma,
    target_index,
    horizon) {
  k <- nrow(sigma)
  ti <- target_index

  ssv <- as_statespace_var(betas, sigma)

  ## Calculate IRFs out to horizon (then adj to 3-dim matrix from DF)
  irf <- irf_ssv(ssv, n_ahead = max(horizon))

  ## Target matrix
  tm <- matrix(0, k, k)
  tm[ti, ti] <- 1

  ## Squared IRF contributions
  contributions <- array(0, dim = c(k, k))
  for (h in min(horizon):max(horizon)) {
    for (i in 1:h) {
      contributions <- contributions + t(irf[, , i]) %*% tm %*% irf[, , i]
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

  return(q)
}
