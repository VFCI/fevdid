#' Impulse Responses.
#'
#' Alias for the function from svars::irf.svars,
#' so that the shock names are correct.
#'
#' @param x SVAR object of class "fevdvar"
#' @param n.ahead Integer specifying the steps.
#' @param ... Currently not used.
#'
#' @return Matrix of forecast error variance decomposition in frequency domain
#' Indexed: frequencies, variables, shocks
#'
#' @rdname irf
#' @name irf
#' @aliases irf.fevdvar
#' @importFrom vars irf
#'
#' @export
#'
#' @examples
#' x <- svars::USA
#' v <- vars::VAR(x, p = 2)
#' mvar <- id_fevdtd(v, "pi", 4:10)
#' vars::irf(mvar)
#'
irf.fevdvar <- function(x, n.ahead = 20, ...) {

    class(x) <- "svars"

    k <- x$K
    variable_names <- colnames(x$y)
    shock_names <- c("Main", paste0("Orth_", 2:k))
    shock_labels <- paste0("epsilon[ ", shock_names, " ]%->%")

    irf <- vars::irf(x, n.ahead = n.ahead, ...)

    names(irf$irf) <- c("V1", apply(
        expand.grid(shock_labels, variable_names), 
        1,
        paste0, collapse = "")
        )

    return(irf)
}
