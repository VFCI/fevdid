test_fe_fev <- function(svar, horizon) {

  fe <- fe(svar, horizon)
  fe_fev_h <- diag(cov(na.omit(fe)))

  fev_fev <- fev(svar, n_ahead = horizon + 1)$fev
  fev_fev_h <- fev_fev[fev_fev$h == horizon, "fev"] |>
    matrix(svar$K, svar$K) |>
    colSums()


  testthat::expect_equal(fe_fev_h, fev_fev_h, tolerance = 0.2)
}

test_that("forecast matches predict function", {

  x <- svars::USA
  x2 <- svars::LN

  var <- vars::VAR(x, p = 2)
  var_const <- vars::VAR(x, p = 2, type = "const")
  var_more_lags <- vars::VAR(x, p = 4, type = "const")
  var_ln <- vars::VAR(x2, p = 4, type = "const")

  svar <- svars::id.chol(var)
  svar_const <- svars::id.chol(var_const)
  svar_more_lags <- svars::id.chol(var_more_lags)
  svar_ln <- svars::id.chol(var_ln)

  for (h in c(1, 2, 3, 5, 10, 20, 50)) {
    test_fe_fev(svar, h)
    test_fe_fev(svar_const, h)
    test_fe_fev(svar_more_lags, h)
    test_fe_fev(svar_ln, h)
  }

})
