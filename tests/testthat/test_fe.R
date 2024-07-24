test_fe_forecast <- function(x, var, horizon) {
  forecast <- forecast(var, horizon = horizon)
  fe_forecast <- rbind(x, matrix(NA, horizon, var$K)) - forecast

  fe <- fe(var, horizon)

  testthat::expect_equal(na.omit(fe_forecast), na.omit(fe), ignore_attr = TRUE, tolerance = 1e-5)
}

test_that("forecast matches predict function", {

  x <- svars::USA
  x2 <- svars::LN

  var <- vars::VAR(x, p = 2)
  var_const <- vars::VAR(x, p = 2, type = "const")
  var_more_lags <- vars::VAR(x, p = 4, type = "const")
  var_ln <- vars::VAR(x2, p = 4, type = "const")

  for (h in c(1, 2, 3, 5, 10, 20, 50)) {
    test_fe_forecast(x, var, h)
    test_fe_forecast(x, var_const, h)
    test_fe_forecast(x, var_more_lags, h)
    test_fe_forecast(x2, var_ln, h)
  }

})
