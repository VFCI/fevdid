#' Identify the main shock by targetting the forecast error
#' variance contribution in the frequency domain.
#'
#' @param var VAR to bootstrap
#' @param id_method function that identifies the Structural VAR
#' @param nboot number of bootstrap iteratons
#' @param horizon length of IRF horizon
#' @param design "recursive" or "fixed".
#' Controls how samples of the data are constructed.
#' @param method "resample" or "wild"
#' @param wild_distr distribution used for "wild" boostrap method.
#' Either "gaussian", "rademacher", or "mammen".
#' @param lower_pctl lower percentile of bootstraps returned in summary.
#' Defaults to 0.16 (for 68% CI).
#' @param upper_pctl upper percentile of bootstraps returned in summary.
#' Defaults to 0.84 (for 68% CI).
#' @param ... Additional arguments passed to the function passed to id_method
#'
#' @return list of bootstrapped VARs and IRFs
#' @export
#'
#' @examples
bootstrap <- function(
    var,
    id_method = NULL,
    nboot = 500,
    horizon = 20,
    design = "recursive",
    method = "resample",
    wild_distr = "gaussian",
    lower_pctl = 0.16,
    upper_pctl = 0.84,
    ...) {
  ## Get VAR objects
  y <- var$y
  p <- var$p
  k <- var$K
  b <- var$B
  a <- var$A_hat
  n <- nrow(y)
  u <- t(stats::resid(var$VAR))


  ## Construct z: stacked y's for statespace var
  z <- t(Reduce(cbind, lapply(1:p, function(l) y[(p + 1 - l):(n - l), ])))

  y_pred <- y[-c(1:p), ] - t(u)

  if (var$type == "const") {
    z <- rbind(1, z)
  } else if (var$type == "none") {
    z <- z
  } else {
    stop("Unsupported VAR type: ", var$type)
  }


  ## Bootstrap the errors
  if (method == "resample") {
    resample <- sample(1:(n - p), size = (n - p) * nboot, replace = TRUE) |>
      matrix(nrow = n - p, ncol = nboot)

    boot_errors <- array(u[, resample], dim = c(k, n - p, nboot))
  } else if (method == "wild") {
    if (wild_distr == "gaussian") {
      adjust <- stats::rnorm(n = (n - p) * nboot)
    } else if (wild_distr == "rademacher") {
      adjust <- sample(c(-1, 1), size = (n - p) * nboot, replace = TRUE)
    } else if (wild_distr == "mammen") {
      cutoff <- (sqrt(5) + 1) / (2 * sqrt(5))
      uniform_draw <- stats::runif(n = (n - p) * nboot, min = 0, max = 1)
      adjust <- ifelse(
        uniform_draw > cutoff, (sqrt(5) + 1) / 2, -(sqrt(5) - 1) / 2
      )
    } else {
      stop("Unsupported wild_distr: ", wild_distr)
    }

    adjust <- matrix(adjust, nrow = n - p, ncol = nboot)
    boot_errors <- array(0, dim = c(k, n - p, nboot))

    for (i in 1:nboot) {
      boot_errors[, , i] <- u * adjust[, i]
    }
  } else {
    stop("Unsupported method.")
  }

  ## Initialize empty arrays to store bootstrap results
  aboots <- array(0, dim = c(dim(a), nboot))
  bboots <- array(0, dim = c(dim(b), nboot))
  irfboots <- list()

  ## Bootsrap Progress Bar
  tick <- 0
  ticks <- 20
  pb <- paste0(c(rep("=", tick), rep("-", ticks - tick)), collapse = "")
  cat("[", pb, "]")

  for (nb in 1:nboot) {
    ## Progrss bar updates
    if (nb %% (nboot / ticks) == 0) {
      tick <- tick + 1
      pb <-
        paste0(c(rep("=", tick), rep("-", ticks - tick)), collapse = "")
      cat("\r[", pb, "]")
    }
    ## Construct data values from the errors
    if (design == "recursive") {
      ystar <- matrix(0, n, k)
      colnames(ystar) <- colnames(y)
      ystar[1:p, ] <- y[1:p, ]

      ## recursively construct next ystar values
      for (i in (p + 1):n) {
        for (j in 1:k) {
          ystar[i, j] <-
            a[j, 1] + a[j, -1] %*% c(t(ystar[(i - 1):(i - p), ])) +
            boot_errors[j, (i - p), nb]
        }
      }

      ystar <- ystar[-c(1:p), ]

      bootvar <- vars::VAR(ystar, p = p, type = var$type)
    } else if (design == "fixed") {
      ystar <- y_pred + t(boot_errors[, , nb])
      bstar <- t(ystar) %*% t(z) %*% solve(z %*% t(z))
      ustar <- ystar - t(bstar %*% z)

      bootvar <- list(
        y = ystar,
        coef_x = bstar,
        residuals = ustar,
        p = p,
        type = var$type
      )
      class(bootvar) <- "var.boot"
    } else {
      stop("Unsupported design: ", design)
    }

    ## Identify bootstrapped VARs
    if (!is.null(id_method)) {
      bootsvar <- id_method(bootvar, ...)
    } else {
      bootsvar <- svars::id.chol(bootvar)
    }

    ## Construct bootstrapped IRFs
    bootirf <- irf(bootsvar, n.ahead = horizon)
    bootirf[[1]]$boot <- nb

    ## Store results
    aboots[, , nb] <- bootsvar$A_hat
    bboots[, , nb] <- bootsvar$B
    irfboots[nb] <- bootirf
  }
  cat("\n")

  ## Calculate summary stats for IRFs
  boot <- h <- name <- variable <- shock <- value <- NULL

  irf_df <- Reduce(rbind, irfboots) |>
    dplyr::rename(h = "V1") |>
    tidyr::pivot_longer(-c(boot, h)) |>
    dplyr::mutate(variable = stringr::str_extract(name, "(?<=%->%).*$")) |>
    dplyr::mutate(shock = stringr::str_extract(name, "(?<= ).*(?= )")) |>
    dplyr::select(!name)

  irf_summarized <- irf_df |>
    dplyr::group_by(h, variable, shock) |>
    dplyr::summarize(
      mean = mean(value),
      median = stats::median(value),
      lower = stats::quantile(value, lower_pctl),
      upper = stats::quantile(value, upper_pctl)
    )

  ## Return bootstrap results
  bootstrap <- list(
    VAR = var,
    A_boots = aboots,
    B_boots = bboots,
    IRF = irf_summarized,
    all_IRF_boots = irfboots,
    id_method = id_method,
    nboot = nboot,
    horizon = horizon,
    design = design,
    method = method,
    wild_distr = wild_distr
  )

  return(bootstrap)
}
