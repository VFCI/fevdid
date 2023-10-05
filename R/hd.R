#' Calculate the historical decomposition for a VAR.
#'
#' @param x A VAR. Currently supports either a 'svars' or 'fevdvar' object.
#' @param ... Not currently used.
#'
#' @return Class 'fevdhd' list of historical decomposition in a data.frame.
#' @export
#'
#' @examples
#' x <- svars::USA
#' v <- vars::VAR(x, p = 2)
#' sv <- svars::id.chol(v)
#' hd <- hd(sv)
#'
hd <- function(
    x,
    ...
    ) {
  UseMethod("hd")
}

#' Method to calculate hd for svars (id.chol)
#'
#' @rdname hd
#' @name hd
#' @aliases hd.svars
#' @export
#'
hd.svars <- function(
    x,
    ...
) {
    k <- x$K
    t <- x$n

    ## Set shock names
    impulse_names <- colnames(x$y)
    response_names <- colnames(x$y)

    ## Set as factors
    impulse_names <-
        factor(impulse_names, levels = impulse_names, ordered = TRUE)
    response_names <-
        factor(response_names, levels = response_names, ordered = TRUE)

    ## Get historical decompositions
    hidec <- lapply(1:k, function(i) {
        svars::hd(x, series = i)$hidec
    })

    ## Pull out just each shock contribution
    hd <- unlist(lapply(hidec, function(i) i[, -(1:2)]))

    ## Pull out the total variation
    total <- unlist(lapply(hidec, function(i) rep(i[, 2], times = k) ))

    ## Tidy to DF
    df <- data.frame(
        t = rep(1:t, times = k * k),
        impulse = rep(impulse_names, each = t, times = k),
        response = rep(response_names, each = t * k),
        hd = hd,
        total = total
    )

    hd <- list(hd = df)
    class(hd) <- "fevdhd"

    return(hd)
}

#' Method to calculate fevdfd for fevdvar (id_fevdfd or id_fevdtd)
#'
#' @param cummulative Boolean.
#' Default to False, set to True for cummulative shocks.
#'
#' @rdname hs
#' @name hs
#' @aliases hs.fevdvar
#' @export
#'
hd.fevdvar <- function(
    x,
    cummulative = FALSE,
    ...
) {
    k <- x$K
    t <- x$n

    ## Set shock names
    impulse_names <- c("Main", paste0("Orth_", 2:k))
    response_names <- colnames(x$y)

    ## Set as factors
    impulse_names <-
        factor(impulse_names, levels = impulse_names, ordered = TRUE)
    response_names <-
        factor(response_names, levels = response_names, ordered = TRUE)

    ## Get historical decompositions
    hidec <- lapply(1:k, function(i) {
        svars::hd(x, series = i)$hidec
    })

    ## Pull out just each shock contribution
    hd <- unlist(lapply(hidec, function(i) i[, -(1:2)]))

    ## Pull out the total variation
    total <- unlist(lapply(hidec, function(i) rep(i[, 2], times = k) ))

    ## Tidy to DF
    df <- data.frame(
        t = rep(1:t, times = k * k),
        impulse = rep(impulse_names, each = t, times = k),
        response = rep(response_names, each = t * k),
        hd = hd,
        total = total
    )

    hd <- list(hd = df)
    class(hd) <- "fevdhd"

    return(hd)
}
