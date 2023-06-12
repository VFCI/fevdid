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
id_fevdtd <- function(var, target, horizon) {

  ## Check parameter values are what is expected
  if (!inherits(var, "varest")) stop("Please pass a VAR from 'vars::VAR'.")

  n <- colnames(var$y)
  k <- length(n)
  ni <- 1:k

  if (is.character(target)){
    if (!target %in% n) stop("Please provide a valid target variable.")
    ti <- which(n %in% target)
    t <- target
  } else if (is.numeric(target)) {
    if (!target %in% ni) stop("Please provide a valid target variable.")
    ti <- target
    t <- n[target]
  }

  if (!is.numeric(horizon)) stop("Please provide an integer valued horizon.")
  if (!all(horizon > 0)) stop("Please provide only positive horizon values.")

  ## Fit a Choleskey SVAR (need orthogonal shocks)
  svar <- svars::id.chol(var)

  ## Calculate IRFs out to horizon
  irf <- svars:::IRF(svar$A_hat[, -1], svar$B, max(horizon)) |>
    simplify2array()

  ## Target matrix
  En <- matrix(0, k, k)
  En[ti, ti] <- 1

  ## Squared IRF contributions
  irf2 <- array(0, dim = c(k, k, max(horizon)))
  for (h in 1:max(horizon)){
    irf2[, , h] <- (max(horizon) + 1 - max(min(horizon),h)) * t(irf[, , h]) %*% En %*% irf[, , h]
  }

  if (max(horizon) == 1) {
    contributions <- irf2[, , 1]
  } else{
    contributions <- rowSums(irf2[, , 1:max(horizon)], dims = 2)
  }


  ## Max eigen value
  e <- eigen(contributions)
  ei <- which(e$value == max(e$value))
  evec <- e$vectors[, ei]

  ## Construct max rotation matrix
  Q <- matrix(0, k, k)
  Q[, 1] <- evec
  Q[, 2:k] <- pracma::nullspace(t(evec))

  ## Insert resulting matrix into var
  mvar <- svar
  mvar$B <- svar$B %*% Q
  mvar$method <- "fevdtd"

  return(mvar)
}
