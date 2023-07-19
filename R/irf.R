#' Impulse Responses.
#'
#' Alias for the function from svars::irf.svars,
#' so that the shock names are correct.
#'
#' @param x SVAR object of class "fevdvar".
#' @param impulse A character vector of the impulses, default is all variables.
#' @param response A character vector of the responses,
#'   default is all variables.
#' @param n.ahead Integer specifying the steps.
#' @param ortho Not used. Here to match generic vars::irf.
#' @param cumulative Not used. Here to match generic vars::irf.
#' @param boot Not used. Here to match generic vars::irf.
#' @param ci Not used. Here to match generic vars::irf.
#' @param runs Not used. Here to match generic vars::irf.
#' @param seed Not used. Here to match generic vars::irf.
#' @param ... Currently not used.
#'
#' @return Matrix of forecast error variance decomposition in frequency domain
#' Indexed: frequencies, variables, shocks
#'
#' @rdname irf
#' @name irf
#' @aliases irf.fevdvar
#' @importFrom  vars irf
#'
#' @export
#'
#' @examples
#' x <- svars::USA
#' v <- vars::VAR(x, p = 2)
#' mvar <- id_fevdtd(v, "pi", 4:10)
#' vars::irf(mvar)
#'
irf.fevdvar <- function(
    x, impulse = NULL, response = NULL, n.ahead = 10,
    ortho = TRUE, cumulative = FALSE, boot = TRUE,
    ci = 0.95, runs = 100, seed = NULL, ...) {
  class(x) <- "svars"

  k <- x$K
  variable_names <- colnames(x$y)
  shock_names <- c("Main", paste0("Orth_", 2:k))
  shock_labels <- paste0("epsilon[ ", shock_names, " ]%->%")

  irf <- vars::irf(x, n.ahead = n.ahead, ...)

  ## Rename the shocks
  names(irf$irf) <- c("V1", apply(
    expand.grid(shock_labels, variable_names),
    1,
    paste0,
    collapse = ""
  ))

  ## Select only the impulses and responses requested
  if (is.null(impulse)) {
    impulse_keep <- rep(TRUE, k)
  } else {
    impulse_keep <- shock_names %in% impulse
  }
  if (is.null(response)) {
    response_keep <- rep(TRUE, k)
  } else {
    response_keep <- variable_names %in% response
  }

  keep <- apply(expand.grid(impulse_keep, response_keep), 1, all)

  irf$irf <- irf$irf[, c(TRUE, keep)]

  return(irf)
}
