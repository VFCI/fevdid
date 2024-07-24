test_forecast_predict <- function(var, horizon) {
  forecast <- forecast(var, horizon = horizon)
  forecast <- forecast[nrow(forecast), ]
  predict <- predict(var, n.ahead = horizon)$fcst |>
    purrr::map(~.x[nrow(.x), "fcst"]) |>
    simplify2array()

  testthat::expect_equal(forecast, predict, ignore_attr = TRUE)
}

test_that("forecast matches predict function", {

  x <- svars::USA
  x2 <- svars::LN

  var <- vars::VAR(x, p = 2)
  var_const <- vars::VAR(x, p = 2, type = "const")
  var_more_lags <- vars::VAR(x, p = 4, type = "const")
  var_ln <- vars::VAR(x2, p = 4, type = "const")

  for (h in c(1, 2, 3, 5, 10, 20, 50)) {
    test_forecast_predict(var, h)
    test_forecast_predict(var_const, h)
    test_forecast_predict(var_more_lags, h)
    test_forecast_predict(var_ln, h)
  }

})
