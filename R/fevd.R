#' Forecast error variance decomposition.
#'
#' Alias for the function from svars::fevd.svars,
#' so that the shock names are correct.
#'
#' @param x SVAR object of class "fevdvar"
#' @param n.ahead Integer specifying the steps.
#' @param ... Currently not used.
#'
#' @return Matrix of forecast error variance decomposition in frequency domain
#' Indexed: frequencies, variables, shocks
#'
#' @rdname fevd
#' @name fevd
#' @aliases fevd.fevdvar
#' @importFrom vars fevd
#'
#' @export
#'
#' @examples
#' x <- svars::USA
#' v <- vars::VAR(x, p = 2)
#' mvar <- id_fevdtd(v, "pi", 4:10)
#' vars::fevd(mvar)
#'
fevd.fevdvar <- function(
  x,
  n.ahead = 10,
  ...
  ) {

  n <- colnames(x$y)
  k <- length(n)

  impulse_names <- n
  if (inherits(x, "fevdvar")) {
    impulse_names <- c("Main", paste0("Orth_", 2:k))
  }
  response_names <- n

  ## Set as factors
  impulse_names <-
    factor(impulse_names, levels = impulse_names, ordered = TRUE)
  response_names <-
    factor(response_names, levels = response_names, ordered = TRUE)

  class(x) <- "svars"
  
  fevd <- vars::fevd(x, n.ahead = n.ahead, ...)

  df <- data.frame(
    h = rep(1:n.ahead, times = k * k),
    impulse = rep(impulse_names, each = n.ahead, times = k),
    response = rep(response_names, each = k * n.ahead),
    fevd = unlist(lapply(fevd, c)),
    row.names = NULL
  )

  fevd <- list(fevd = df)
  class(fevd) <- c("fevdfevd")

  return(fevd)
}
