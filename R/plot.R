#' Plot the forecast variance decomposition.
#'
#' Alias for the function plot in frequency domain.
#'
#' @param x object of class "fevdfd"
#' @param y Not used.
#' @param stacked Boolean.
#' True for stacked columns, False (default) for unstacked line chart.
#' @param vlines Vector of x-axis points at which to draw a vline.
#' Useful for highlighting areas.
#' @param ... Currently not used.
#'
#' @return ggplot of fevdfd
#'
#' @rdname plot
#' @name plot
#' @aliases plot.fevdfd
#'
#' @export
#'
plot.fevdfd <- function(x, y, stacked = FALSE, vlines = NULL, ...) {

    ## Declare data.frame variables so check() doesn't complain
    f <- shocks <- variable <- value <- NULL

    plot_data <- Reduce(rbind, lapply(names(x), function(n) {
        data <- unlist(x[[n]][, -1], use.names = FALSE)
        f <- x[[n]]$f
        shocks <- names(x[[n]][, -1])

        data.frame(
            f = rep(f, length(shocks)),
            shocks = rep(shocks, each = length(f)),
            variable = n,
            value = data
        )
    }))

    plot <- plot_data |>
        ggplot2::ggplot(ggplot2::aes(
            x = f,
            y = value,
            color = shocks,
            fill = shocks
        )) +
        ggplot2::facet_wrap(
            ggplot2::vars(variable), ncol = 1, scales = "free_y") +
            ggplot2::theme_bw()

    if (stacked == FALSE) {
        plot <- plot +
            ggplot2::geom_line()
    } else {
        plot <- plot +
            ggplot2::geom_col(position = "stack")
    }

    if (!is.null(vlines)) {
        plot <- plot +
            ggplot2::geom_vline(xintercept = vlines)
    }

    return(plot)
}

#'
#' Plot the forecast variance decomposition in time domain.
#'
#' @param x object of class "svarfevd"
#' @param y Not used.
#' @param stacked Boolean.
#' True for stacked columns, False (default) for unstacked line chart.
#' @param vlines Vector of x-axis points at which to draw a vline.
#' Useful for highlighting areas.
#' @param ... Currently not used.
#'
#' @return ggplot of fevdfd
#'
#' @rdname plot
#' @name plot
#' @aliases plot.svarfevd
#'
#' @export
#'
plot.svarfevd <- function(x, y, stacked = FALSE, vlines = NULL, ...) {

    ## Declare data.frame variables so check() doesn't complain
    h <- shocks <- variable <- value <- NULL

    plot_data <- Reduce(rbind, lapply(names(x), function(n) {
        data <- unlist(x[[n]], use.names = FALSE)
        h <- seq_len(nrow(x[[n]]))
        shocks <- names(x[[n]])

        data.frame(
            h = rep(h, length(shocks)),
            shocks = rep(shocks, each = length(h)),
            variable = n,
            value = data
        )
    }))

    plot <- plot_data |>
        ggplot2::ggplot(ggplot2::aes(
            x = h,
            y = value,
            color = shocks,
            fill = shocks
        )) +
        ggplot2::facet_wrap(
            ggplot2::vars(variable), ncol = 1, scales = "free_y") +
            ggplot2::theme_bw()

    if (stacked == FALSE) {
        plot <- plot +
            ggplot2::geom_line()
    } else {
        plot <- plot +
            ggplot2::geom_col(position = "stack")
    }

    if (!is.null(vlines)) {
        plot <- plot +
            ggplot2::geom_vline(xintercept = vlines)
    }

    return(plot)
}
