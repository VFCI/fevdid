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

  ## Calculate IRFs out to horizon (then adj to 3-dim matrix from DF)
  irf <- vars::irf(svar, n.ahead = max(horizon))[[1]][, -1] |>
    apply(1, matrix, simplify = FALSE, nrow = k, ncol = k, byrow = TRUE) |>
    simplify2array()

  ## Target matrix
  tm <- matrix(0, k, k)
  tm[ti, ti] <- 1

  ## Squared IRF contributions
  irf2 <- array(0, dim = c(k, k, max(horizon)))
  for (h in 1:(max(horizon))) {
    h_weight <- max(horizon) + 1 - max(min(horizon), h)
    irf2[, , h] <- h_weight * t(irf[, , h]) %*% tm %*% irf[, , h]
  }

  if (max(horizon) == 1) {
    contributions <- irf2[, , 1]
  } else {
    contributions <- rowSums(irf2[, , 1:max(horizon)], dims = 2)
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
  mvar$method <- "fevdtd"

  return(mvar)
}
