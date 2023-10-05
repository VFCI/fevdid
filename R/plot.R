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
plot.fevdfd <- function(x, y, stacked = TRUE, vlines = NULL, ...) {
  ## Declare data.frame variables so check() doesn't complain
  f <- impulse <- response <- fevdfd <- NULL

  plot_data <- x$fevdfd

  plot <- plot_data |>
    ggplot2::ggplot(ggplot2::aes(
      x = f,
      y = fevdfd,
      color = impulse,
      fill = impulse
    )) +
    ggplot2::facet_wrap(
      ggplot2::vars(response),
      ncol = 1, scales = "free_y"
    ) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_hue() +
    ggplot2::scale_fill_hue() +
    ggplot2::scale_x_continuous(limits = c(0, pi))

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


#' Plot the forecast variance decomposition.
#'
#' Alias for the function plot in frequency domain.
#'
#' @param x object of class "fevfd"
#' @param y Not used.
#' @param stacked Boolean.
#' True for stacked columns, False (default) for unstacked line chart.
#' @param vlines Vector of x-axis points at which to draw a vline.
#' Useful for highlighting areas.
#' @param ... Currently not used.
#'
#' @return ggplot of fevfd
#'
#' @rdname plot
#' @name plot
#' @aliases plot.fevfd
#'
#' @export
#'
plot.fevfd <- function(x, y, stacked = TRUE, vlines = NULL, ...) {
  ## Declare data.frame variables so check() doesn't complain
  f <- impulse <- response <- fevfd <- NULL

  plot_data <- x$fevfd

  plot <- plot_data |>
    ggplot2::ggplot(ggplot2::aes(
      x = f,
      y = fevfd,
      color = impulse,
      fill = impulse
    )) +
    ggplot2::facet_wrap(
      ggplot2::vars(response),
      ncol = 1, scales = "free_y"
    ) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_hue() +
    ggplot2::scale_fill_hue() +
    ggplot2::scale_x_continuous(limits = c(0, pi))

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
#' @param x object of class "fevdfevd"
#' @param y Not used.
#' @param stacked Boolean.
#' True for stacked columns, False (default) for unstacked line chart.
#' @param vlines Vector of x-axis points at which to draw a vline.
#' Useful for highlighting areas.
#' @param ... Currently not used.
#'
#' @return ggplot of fevd
#'
#' @rdname plot
#' @name plot
#' @aliases plot.fevdfevd
#'
#' @export
#'
plot.fevdfevd <- function(x, y, stacked = TRUE, vlines = NULL, ...) {
  ## Declare data.frame variables so check() doesn't complain
  h <- impulse <- response <- value <- NULL

  plot_data <- x$fevd

  plot <- plot_data |>
    ggplot2::ggplot(ggplot2::aes(
      x = h,
      y = fevd,
      color = impulse,
      fill = impulse
    )) +
    ggplot2::facet_wrap(
      ggplot2::vars(response),
      ncol = 1, scales = "free_y"
    ) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_hue() +
    ggplot2::scale_fill_hue()

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
#' Plot the Impulse Response Functions
#'
#' @param x object of class "fevdirf"
#' @param y Not used.
#' @param impulse_as Default to "colors", with all shocks on the same facet.
#' Other option is "cols" which puts each shock as its own column facet.
#' @param ... Currently not used.
#'
#' @return ggplot of irf
#'
#' @rdname plot
#' @name plot
#' @aliases plot.fevdirf
#'
#' @export
#'
plot.fevdirf <- function(x, y = NULL, impulse_as = "colors", ...) {
  ## Declare data.frame variables so check() doesn't complain
  h <- irf <- impulse <- response <- NULL

  plot_data <- x$irf

  plot <- plot_data |>
    ggplot2::ggplot(ggplot2::aes(
      x = h,
      y = irf
    )) +
    ggplot2::geom_hline(yintercept = 0, size = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_hue() +
    ggplot2::scale_fill_hue()

  if (impulse_as == "colors") {
    plot <- plot +
      ggplot2::geom_line(
        ggplot2::aes(color = impulse)
      ) +
      ggplot2::facet_wrap(
        ggplot2::vars(response),
        ncol = 1, scales = "free_y"
      )
  } else if (impulse_as == "cols") {
    plot <- plot +
      ggplot2::geom_line() +
      ggplot2::facet_grid(
        rows = ggplot2::vars(response),
        cols = ggplot2::vars(impulse),
        scales = "free_y"
      )
  } else {
    stop("Please set imuplse to 'colors' or 'cols'.")
  }

  return(plot)
}


#'
#' Plot the forecast error variance in the time domain
#'
#' @param x object of class "fevdfev"
#' @param y Not used.
#' @param impulse_as Default to "colors", with all shocks on the same facet.
#' Other option is "cols" which puts each shock as its own column facet.
#' @param ... Currently not used.
#'
#' @return ggplot of fevdfev
#'
#' @rdname plot
#' @name plot
#' @aliases plot.fevdfev
#'
#' @export
#'
plot.fevdfev <- function(x, y = NULL, impulse_as = "colors", ...) {
  ## Declare data.frame variables so check() doesn't complain
  h <- fev <- impulse <- response <- NULL

  plot_data <- x$fev

  plot <- plot_data |>
    ggplot2::ggplot(ggplot2::aes(
      x = h,
      y = fev
    )) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_hue() +
    ggplot2::scale_fill_hue()

  if (impulse_as == "colors") {
    plot <- plot +
      ggplot2::geom_line(
        ggplot2::aes(color = impulse)
      ) +
      ggplot2::facet_wrap(
        ggplot2::vars(response),
        ncol = 1, scales = "free_y"
      )
  } else if (impulse_as == "cols") {
    plot <- plot +
      ggplot2::geom_line() +
      ggplot2::facet_grid(
        rows = ggplot2::vars(response),
        cols = ggplot2::vars(impulse)
      )
  } else {
    stop("Please set imuplse to 'colors' or 'cols'.")
  }

  return(plot)
}

#'
#' Plot the historical shocks.
#'
#' @param x object of class "fevdhs"
#' @param y Not used.
#' @param impulse_as Default to "colors", with all shocks on the same facet.
#' Other option is "cols" which puts each shock as its own column facet.
#' @param ... Currently not used.
#'
#' @return ggplot of hs
#'
#' @rdname plot
#' @name plot
#' @aliases plot.fevdhs
#'
#' @export
#'
plot.fevdhs <- function(x, y = NULL, impulse_as = "colors", ...) {
  ## Declare data.frame variables so check() doesn't complain
  t <- hs <- impulse <- NULL

  plot_data <- x$hs

  plot <- plot_data |>
    ggplot2::ggplot(ggplot2::aes(
      x = t,
      y = hs
    )) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_hue() +
    ggplot2::scale_fill_hue()

    if (impulse_as == "colors") {
      plot <- plot +
      ggplot2::geom_line(
        ggplot2::aes(color = impulse)
      )
    } else if (impulse_as == "cols") {
      plot <- plot +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(
        ggplot2::vars(impulse),
        ncol = 1
      )
    } else {
      stop("Please set imuplse to 'colors' or 'cols'.")
    }

  return(plot)
}

#'
#' Plot the historical decomposition.
#'
#' @param x object of class "fevdhd"
#' @param y Not used.
#' @param ... Currently not used.
#'
#' @return ggplot of hd
#'
#' @rdname plot
#' @name plot
#' @aliases plot.fevdhd
#'
#' @export
#'
plot.fevdhd <- function(x, y = NULL, ...) {
  ## Declare data.frame variables so check() doesn't complain
  t <- hd <- impulse <- response <- total <- NULL

  plot_data <- x$hd

  plot <- plot_data |>
    ggplot2::ggplot(ggplot2::aes(
      x = t,
    )) +
    ggplot2::geom_col(ggplot2::aes(
      y = hd,
      fill = impulse
    )) +
    ggplot2::geom_line(ggplot2::aes(
      y = total,
      color = "Total"
    )) +
    ggplot2::facet_wrap(
      ggplot2::vars(response),
      ncol = 1
    ) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_hue() +
    ggplot2::scale_color_manual(
      values = c(Total = "black")
      ) +
    ggplot2::labs(color = NULL)

  return(plot)
}