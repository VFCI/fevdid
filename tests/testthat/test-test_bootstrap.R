
test_wildboot <- function(p, distr, design, nboot, tols) {
  x <- svars::USA
  v <- vars::VAR(x, p = p, type = "const")

  sv <- svars::id.chol(v)

  bv <- svars::wild.boot(
    sv, design = design, distr = distr, n.ahead = 10, nboot = nboot)

  bv2 <- bootstrap(sv, svars::id.chol, n_ahead = 10, nboot = nboot,
    design = design, method = "wild", wild_distr = distr
    )

  compare_boots(bv, bv2, tols)
}

compare_boots <- function(bv, bv2, tols) {
  expect_equal(bv$A_hat_boot_mean, bv2$A_mean, tolerance = tols[1])

  bv_boots <- lapply(bv$bootstrap, function(x) {
    x[[1]][, -1] |>
    t() |>
    unlist() |>
    array(dim = c(3, 3, 10)) |>
    aperm(c(2, 1, 3))
    })
  bv_boots <- array(unlist(bv_boots), dim = c(3, 3, 10, bv$nboot))
  bv_mean_boot <- rowMeans(bv_boots, dims = 3)
  bv2_mean_boot <- bv2$IRF_mean

  bv_mean_boot <- array(aperm(bv_mean_boot, c(3, 2, 1)), dim = c(10, 9))
  bv2_mean_boot <- array(aperm(bv2_mean_boot, c(3, 2, 1)), dim = c(10, 9))

  ## Compare equality IRF by IRF
  expect_equal(
    bv_mean_boot,
    bv2_mean_boot,
    tolerance = tols[2]
    )
}


test_that("bootstrap matches svars::wild.boot", {

  ## Test Gaussian Wild Distribution with Recursive
  test_wildboot(2, "gaussian", "recursive", 250, c(0.1, 0.1))

  ## Test Rademacher Wild Distribution with Recursive
  test_wildboot(2, "rademacher", "recursive", 250, c(0.1, 0.1))

  ## Test mammen Wild Distribution with Recursive
  test_wildboot(2, "mammen", "recursive", 250, c(0.1, 0.1))

  ## Test Gaussian Wild Distribution with Fixed
  test_wildboot(2, "gaussian", "fixed", 250, c(0.1, 0.1))

  ## Test different lags
  test_wildboot(1, "gaussian", "recursive", 250, c(0.1, 0.1))
  test_wildboot(5, "gaussian", "recursive", 250, c(0.1, 0.1))
  test_wildboot(10, "gaussian", "recursive", 250, c(0.1, 0.1))
})


test_that("bootstrap with resample is close to svars::wild.boot", {
  x <- svars::USA
  v <- vars::VAR(x, p = 2, type = "const")

  sv <- svars::id.chol(v)

  bv <- svars::wild.boot(
    sv, design = "recursive", distr = "gaussian", n.ahead = 10, nboot = 250
    )

  bv2 <- bootstrap(
    sv, svars::id.chol, n_ahead = 10, nboot = 250,
    design = "recursive", method = "resample", wild_distr = "gaussian"
    )

  ## Bigger tolerance because these are expected to be kind of different
  compare_boots(bv, bv2, c(0.25, 0.1))

})
