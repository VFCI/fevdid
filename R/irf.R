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
#' @param as_matrix Default to False.
#' Set to true to return a 3D matrix instead of a data.frame.
#' @param ... Currently not used.
#'
#' @return Matrix of forecast error variance decomposition in frequency domain
#' Indexed: frequencies, variables, shocks
#'
#' @rdname irf
#' @name irf
#' @aliases irf.fevdvar
#' @importFrom  vars irf
#' @importFrom  expm %^%
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
    x,
    impulse = NULL,
    response = NULL,
    n.ahead = 10,
    ortho = TRUE,
    cumulative = FALSE,
    boot = TRUE,
    ci = 0.95,
    runs = 100,
    seed = NULL,
    as_matrix = FALSE,
    ...) {
  class(x) <- "svars"

  k <- x$K
  response_names <- colnames(x$y)
  impulse_names <- x$impulse_names


  ## Represent as state space var
  ssv <- as_statespace_var(x$A_hat, x$B)
  irf <- irf_ssv(ssv, n_ahead = n.ahead)

  if (as_matrix) {
    return(irf)
  }

  ## Create data.frame
  irf_df <- data.frame(
    h = rep(1:n.ahead, each = k * k),
    impulse = rep(impulse_names, each = k, times = n.ahead),
    response = rep(response_names, times = k * n.ahead),
    irf = c(irf)
  )

  ## Select only the impulses and responses requested
  if (!is.null(impulse)) {
    irf_df <- irf_df[irf_df$impulse %in% impulse, ]
  }
  if (!is.null(response)) {
    irf_df <- irf_df[irf_df$response %in% response, ]
  }

  irf <- list(irf = irf_df)
  class(irf) <- "fevdirf"


  return(irf)
}


#' @rdname irf
#' @name irf
#' @aliases irf.bvartools
#'
#' @export
#'
irf.bvartools <- function(
    x,
    impulse = NULL,
    response = NULL,
    n.ahead = 10,
    ortho = TRUE,
    cumulative = FALSE,
    boot = TRUE,
    ci = 0.95,
    runs = 100,
    seed = NULL,
    ...) {
  k <- nrow(x$y)
  var_names <- rownames(x$y)
  iterations <- ncol(x$Sigma)

  if (x$method %in% c("id_fevdfd", "id_fevdtd")) {
    shock_names <- c("Main", paste0("Orth_", 2:k))
  } else {
    shock_names <- paste0("V", 1:k)
  }

  iter_irfs <- array(NA, dim = c(k, k, n.ahead, iterations))

  for (i in 1:iterations) {
    a <- matrix(x$A[, i], k)
    c <- matrix(x$C[, i], k, 1)

    a_hat <- cbind(c, a)

    sig <- matrix(x$Sigma[, i], k, k)

    ssv <- as_statespace_var(a_hat, sig)

    iter_irfs[, , , i] <- irf_ssv(ssv, n_ahead = n.ahead)
  }

  ## Create tidy IRF DF
  irf_df <- data.frame(
    h = rep(1:n.ahead, each = k * k),
    impulse = rep(shock_names, each = k, times = n.ahead),
    response = rep(var_names, times = k * n.ahead),
    mean = c(rowMeans(iter_irfs, dims = 3)),
    median = c(apply(iter_irfs, c(1, 2, 3), stats::median)),
    lower = c(apply(iter_irfs, c(1, 2, 3), stats::quantile, probs = 0.16)),
    upper = c(apply(iter_irfs, c(1, 2, 3), stats::quantile, probs = 0.84))
  )

  return(irf_df)
}


#'
#' Simple IRF function for statespace VARs
#'
#' @param ssv Statespace VAR class.
#' @param n_ahead integer. How far out to calculate the impulse response.
#'
irf_ssv <- function(
    ssv,
    n_ahead = 10) {
  k <- nrow(ssv$my)

  irf <- array(0, dim = c(k, k, n_ahead))

  for (h in seq_len(n_ahead)) {
    irf[, , h] <- ssv$my %*% (ssv$mx %^% (h - 1) %*% ssv$me)
  }

  return(irf)
}
