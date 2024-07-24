#' Calculate the forecasts for a VAR out to a specified horizon
#' Conditions on information at time "t"
#'
#' @param var vars::VAR or svars object
#' @param horizon number of steps out to calculate forecast
#'
#' @return forecast
#'
#' @export
#'
#' @importFrom  expm %^%
#'
#' @examples
#' x <- svars::USA
#' v <- vars::VAR(x, p = 2)
#' forecast(v, 10)
#'
forecast <- function(var, horizon = 1) {
  ## Declare this to avoid "no visible binding" warning
  `%^%` <- expm::`%^%`

  if (inherits(var, "varest")) {
    svar <- svars::id.chol(var)
  } else if (inherits(var, "svars")) {
    svar <- var
  }

  ssv <- as_statespace_var(svar$A_hat, svar$B)

  x <- var$datamat[, 1:var$K]

  x_lag1 <- var$datamat[1, (var$K + 1):(var$K * (var$p + 1))]
  x_lag2 <- var$datamat[, 1:(var$K * var$p)]
  names(x_lag1) <- names(x_lag2)

  x_lag <- rbind(x_lag1, x_lag2)

  if (var$type == "const") {
    const_a <-
      purrr::map(0:(horizon - 1), ~ ssv$mx %^% .x) |>
      simplify2array() |>
      rowSums(dims = 2)

    const <- const_a %*% c(svar$A_hat[, "const"], rep(0, var$K * (var$p - 1)))
  } else {
    const <- matrix(0, var$K * var$p, 1)
  }

  non_const <- (ssv$mx %^% (horizon)) %*% t(as.matrix(x_lag))

  x_hat_h <- t(c(ssv$my %*% const) + ssv$my %*% non_const)

  x_hat_h <- rbind(matrix(NA, var$p + (horizon - 1), var$K), x_hat_h)
  colnames(x_hat_h) <- colnames(x)

  return(x_hat_h)
}
