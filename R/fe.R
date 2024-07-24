#' Calculate the historical forecast error for a VAR out to a specified horizon
#'
#' @param var vars::VAR or svars object
#' @param horizon number of steps out to calculate fe
#'
#' @return forecast error
#'
#' @export
#'
#' @importFrom  expm %^%
#'
#' @examples
#' x <- svars::USA
#' v <- vars::VAR(x, p = 2)
#' fe(v, 10)
#'
fe <- function(var, horizon = 1) {

  ## Declare this to avoid "no visible binding" warning
  `%^%` <- expm::`%^%`

  if (inherits(var, "varest")) {
    svar <- svars::id.chol(var)
    res <- stats::residuals(var)
  } else if (inherits(var, "svars")) {
    svar <- var
    res <- stats::residuals(svar$VAR)
  }

  ssv <- as_statespace_var(svar$A_hat, svar$B)
  res <- cbind(res, matrix(0, nrow(res), var$K * (var$p - 1)))

  fe <- matrix(NA, nrow(res), var$K)
  names(fe) <- names(res)

  for (t in seq_len(nrow(res) - (horizon - 1))) {
    fe_t <- 0

    for (h in seq_len(horizon)) {
      irf_h <- ssv$mx %^% (h - 1)
      res_h <- res[t + (horizon - 1) - (h - 1), ]

      e <- irf_h %*% res_h

      fe_t <- fe_t + e
    }
    fe[t + (horizon - 1), ] <- ssv$my %*% fe_t
  }

  return(fe)
}
