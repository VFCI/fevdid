#' Calculate the forecast error variance for a VAR out to a specified horizon
#'
#' @param var vars::VAR or svars object
#' @param h number of steps out to calculate fev
#'
#' @return forecast error variance
#'
fev <- function(var, h) {
  n <- colnames(var$y)
  k <- length(n)

  ## Calculate IRFs out to horizon (then adj to 3-dim matrix from DF)
  if (inherits(var, "fevdvar")) {
    irf <- vars::irf(var, n.ahead = h, as_matrix = TRUE)
  } else {
    irf <- vars::irf(var, n.ahead = h)[[1]][, -1] |>
      apply(1, matrix, simplify = FALSE, nrow = k, ncol = k, byrow = TRUE) |>
      simplify2array()
  }


  # Forecast Error Variance calculation
  fe <- list()
  for (i in 1:k) {
    fe[[i]] <- as.data.frame(t(irf[i, , ]))
    colnames(fe[[i]]) <- n
  }
  names(fe) <- n
  fe2 <- fe

  for (i in seq_along(fe)) {
    for (j in 1:h) {
      fe2[[i]][j, ] <- (colSums(fe[[i]][j:1, ]^2))
    }
  }

  return(fe2)
}
