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
fevd.fevdvar <- function(x, n.ahead = 10, ...) {
  class(x) <- "svars"

  k <- x$K
  fevd <- vars::fevd(x, n.ahead = n.ahead, ...)

  for (i in seq_along(fevd)) {
    colnames(fevd[[i]]) <- c("Main", paste0("Orth_", 2:k))
  }

  class(fevd) <- c("fevdvarfevd", "svarfevd")

  return(fevd)
}
