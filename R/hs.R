#' Calculate the historical structural shocks for a VAR.
#'
#' @param x A VAR. Currently supports either a 'svars' or 'fevdvar' object.
#' @param ... Not currently used.
#'
#' @return Class 'fevdhs' list of historical shocks in a data.frame.
#' @export
#'
#' @examples
#' x <- svars::USA
#' v <- vars::VAR(x, p = 2)
#' sv <- svars::id.chol(v)
#' hs <- hs(sv)
#'
hs <- function(
    x,
    ...) {
  UseMethod("hs")
}

#' Method to calculate hs for svars (id.chol)
#'
#' @param cumulative Boolean.
#' Default to False, set to True for cummulative shocks.
#'
#' @rdname hs
#' @name hs
#' @aliases hs.svars
#' @export
#'
hs.svars <- function(
    x,
    cumulative = FALSE,
    ...) {
  k <- x$K

  ## Set shock names
  impulse_names <- colnames(x$y)

  ## Set as factors
  impulse_names <-
    factor(impulse_names, levels = impulse_names, ordered = TRUE)

  ## Get residuals and Sigma
  residuals <- stats::resid(x$VAR)
  sigma <- x$B
  t <- nrow(residuals)

  hs <- hist_shocks(residuals, sigma)

  ## Tidy to DF
  df <- data.frame(
    t = rep(1:t, times = k),
    impulse = rep(impulse_names, each = t),
    hs = c(hs)
  )

  ## Cummulative shocks
  if (cumulative == TRUE) {
    df <- data.frame(
      t = rep(1:t, times = k),
      impulse = rep(impulse_names, each = t),
      hs = c(apply(hs, 2, cumsum))
    )
  }

  hs <- list(hs = df)
  class(hs) <- "fevdhs"

  return(hs)
}

#' Method to calculate hs for fevdvar (id_fevdfd or id_fevdtd)
#'
#' @param cumulative Boolean.
#' Default to False, set to True for cummulative shocks.
#'
#' @rdname hs
#' @name hs
#' @aliases hs.fevdvar
#' @export
#'
hs.fevdvar <- function(
    x,
    cumulative = FALSE,
    ...) {
  k <- x$K

  ## Set shock names
  impulse_names <- c("Main", paste0("Orth_", 2:k))

  ## Set as factors
  impulse_names <-
    factor(impulse_names, levels = impulse_names, ordered = TRUE)

  ## Get residuals and Sigma
  residuals <- stats::resid(x$VAR)
  sigma <- x$B
  t <- nrow(residuals)

  hs <- hist_shocks(residuals, sigma)

  ## Tidy to DF
  df <- data.frame(
    t = rep(1:t, times = k),
    impulse = rep(impulse_names, each = t),
    hs = c(hs)
  )

  ## Cummulative shocks
  if (cumulative == TRUE) {
    df <- data.frame(
      t = rep(1:t, times = k),
      impulse = rep(impulse_names, each = t),
      hs = c(apply(hs, 2, cumsum))
    )
  }

  hs <- list(hs = df)
  class(hs) <- "fevdhs"

  return(hs)
}


#' Find the historical shocks for a VAR
#'
#' @param residuals Matrix of empirical VAR residuals
#' @param sigma Matrix mapping empirical residuals to structural
#'
#' @return historical shocks, matrix
#'
hist_shocks <- function(
    residuals,
    sigma) {
  hs <- residuals %*% sigma

  return(hs)
}
