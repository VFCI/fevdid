#' Identify a VAR with a Cholesky ordering assumption.
#' Assumes the given order of the variables in the VAR,
#' reorder the data first to get a different ordering.
#'
#' @param x a vars::VAR
#'
#' @return SVAR with cholesky id
#' @export
#'
id_ordered_chol <- function(
  x
) {

  n <- colnames(x$y)
  k <- length(n)

  chol_var <- svars::id.chol(x)

  ## svars package does a degree of freedom adjustment, so recalc the B matrix
  b <- t(chol(stats::cov(stats::residuals(x))))

  ## q is the "rotation" matrix, in this case, just the identity
  q <- diag(x$K)

  chol_var$Q <- q
  chol_var$B <- b %*% q
  chol_var$method <- "id_ordered_chol"
  chol_var$impulse_names <- paste0("Chol_", 1:k)

  class(chol_var) <- c("fevdvar", "svars")

  return(chol_var)
}