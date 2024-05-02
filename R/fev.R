#' Calculate the forecast error variance for a VAR out to a specified horizon
#'
#' @param var vars::VAR or svars object
#' @param n_ahead number of steps out to calculate fev
#'
#' @return forecast error variance
#'
#' @export
#'
#' @examples
#' x <- svars::USA
#' v <- vars::VAR(x, p = 2)
#' mvar <- id_fevdtd(v, "pi", 4:10)
#' fev(mvar)
#'
fev <- function(var, n_ahead = 20) {
  n <- colnames(var$y)
  k <- length(n)

  impulse_names <- n
  if (inherits(var, "fevdvar")) {
    impulse_names <- var$impulse_names
  }
  response_names <- n

  ## Set as factors
  impulse_names <-
    factor(impulse_names, levels = impulse_names, ordered = TRUE)
  response_names <-
    factor(response_names, levels = response_names, ordered = TRUE)

  ## Calculate IRFs out to horizon (then adj to 3-dim matrix from DF)
  if (inherits(var, "fevdvar")) {
    irf <- vars::irf(var, n.ahead = n_ahead, as_matrix = TRUE)
  } else {
    irf <- vars::irf(var, n.ahead = n_ahead)[[1]][, -1] |>
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
    for (j in 1:n_ahead) {
      fe2[[i]][j, ] <- (colSums(fe[[i]][j:1, ]^2))
    }
  }

  ## Tidy to a DF
  df <- data.frame(
    h = rep(1:n_ahead, times = k * k),
    impulse = rep(impulse_names, each = n_ahead, times = k),
    response = rep(response_names, each = k * n_ahead),
    fev = unlist(lapply(fe2, c)),
    row.names = NULL
  )

  fev <- list(fev = df)
  class(fev) <- "fevdfev"

  return(fev)
}
