test_that("Check that fevdtd > id.chol shocks", {

  x <- svars::USA

  v <- vars::VAR(x, p = 2)
  svar <- svars::id.chol(v)
  svar$B <- t(chol(cov(residuals(v))))

  h <- list(1, 1:10, 1:20, 1:40, 1:60, 5:10, 10:15, 15:20, 20:40, 40:60)
  t <- c("x", "pi", "i")
  iter <- expand.grid(t = t, h = h)

  for (i in seq_len(nrow(iter))) {
    ti <- iter[[i, "t"]]
    hi <- iter[[i, "h"]]

    mvar <- id_fevdtd(v, ti, hi)

    fevdm <- fev(mvar, h = max(hi))
    fevds <- fev(svar, h = max(hi))

    m <- sum(fevdm[[ti]][, 1][hi])
    m2 <- colSums(fevds[[ti]][hi, ])

    expect_true(
      all(m >= m2),
      label = paste0("Iter: ", i, "\n", m, "\n", paste(m2, collapse = ", "))
      )
  }

})


test_that("Check parameter validation", {

  x <- svars::USA
  v <- vars::VAR(x, p = 2)

  expect_error(id_fevdtd(v, "AAAA", 1))
  expect_error(id_fevdtd(v, 0, 1))
  expect_error(id_fevdtd(v, "x", 0))
  expect_error(id_fevdtd(v, "x", "pi"))
  expect_error(id_fevdtd(v, 4, 1))

})
