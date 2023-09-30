#' Convert a VAR with p lags to a statespace representation with just one lag
#'
#' @param betas matrix of coefficients from a VAR.
#' Should have constant coefficient in the first column if included.
#' Dimensions should be k x (1 + k * p), or k x (k x p) if no constant.
#' @param sigma matrix of the var-cov of the residuals.
#' Should be k x k in dimension.
#'
#' @return matrix coefficients for VAR 1 representation
#'

as_statespace_var <- function(betas, sigma) {
  k <- nrow(sigma)
  p <- floor(ncol(betas) / k)

  ## Get coeefficient matrix, ignore constant
  if ((p * k) + 1 == ncol(betas)) {
    a <- betas[, -1]
  } else if (p * k == ncol(betas)) {
    a <- betas
  } else {
    stop("Betas have the wrong dimensions.")
  }

  ## Create Var(1) objects
  my <- cbind(diag(k), matrix(0, k, k * (p - 1)))
  mx <- rbind(a, cbind(diag(k * (p - 1)), matrix(0, k * (p - 1), k)))
  me <- rbind(sigma, matrix(0, k * (p - 1), k))

  if (p == 1) my <- diag(k)

  ssv <- list(
    my = my,
    mx = mx,
    me = me
  )

  class(ssv) <- "statespacevar"

  return(ssv)
}
