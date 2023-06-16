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
  svar$B <- t(chol(cov(residuals(var))))

  ## Calculate IRFs out to horizon (then adj to 3-dim matrix from DF)
  irf <- vars::irf(svar, n.ahead = max(horizon))[[1]][, -1] |>
    apply(1, matrix, simplify = FALSE, nrow = k, ncol = k, byrow = TRUE) |>
    simplify2array()

  vtmp <- t(irf[ti, , ])

  V0 <- matrix(0, k, k)
  V1 <- 0
  for (i in 1:max(horizon)){
    V0 <- V0 + vtmp[i, ] %*% t(vtmp[i, ])
    V1 <- V1 + t(vtmp[i, ]) %*% vtmp[i, ]
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
