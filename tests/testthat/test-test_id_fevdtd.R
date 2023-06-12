test_that("Check that fevdtd > id.chol shocks", {
  h <- 4:10
  x <- svars::USA
  v <- vars::VAR(x, p = 2)
  svar <- svars::id.chol(v)
  mvar <- id_fevdtd(v, "pi", h)
  d1 <- svars:::fevd.svars(mvar, n.ahead = max(h))
  d2 <- svars:::fevd.svars(svar, n.ahead = max(h))
  m <- mean(d1[[2]][, 1][h])
  m2 <- colMeans(d2[[2]][h, ])
  expect_true(all(m > m2))
})
