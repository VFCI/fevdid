#' Identify the main shock by targetting the forecast error
#' variance contribution in the frequency domain.
#'
#' @param var VAR to bootstrap
#' @param id_method function that identifies the Structural VAR
#' @param nboot number of bootstrap iteratons
#' @param n_ahead length of IRF horizon
#' @param design "recursive" or "fixed".
#' Controls how samples of the data are constructed.
#' @param method "resample" or "wild"
#' @param wild_distr distribution used for "wild" boostrap method.
#' Either "gaussian", "rademacher", or "mammen".
#' @param bias_adjust Boolean.  Set to TRUE to calculate the bias
#' and run the bootstrap again.
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
    n_ahead = 20,
    design = "recursive",
    method = "resample",
    wild_distr = "gaussian",
    bias_adjust = FALSE,
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

  ## Variable and Shock names
  variable_names <- colnames(y)
  if (inherits(var, "fevdvar")){
    shock_names <-  c("Main", paste0("Orth_", 2:k))
  } else {
    shock_names <- colnames(y)
  }

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
  irfboots <- array(0, dim = c(k, k, n_ahead, nboot))

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
      bootsvar <- bootvar
    }

    ## Construct bootstrapped IRFs
    bootirf_df <- irf(bootsvar, n.ahead = n_ahead)[[1]]
    bootirf <- array(unlist(t(bootirf_df[, -1])), dim = c(k, k, n_ahead))

    ## Store results
    aboots[, , nb] <- bootsvar$A_hat
    bboots[, , nb] <- bootsvar$B
    irfboots[, , , nb] <- bootirf
  }
  cat("\n")

  ## Calculate summary stats for IRFs
  a_mean <- rowMeans(aboots, dims = 2)
  b_mean <- rowMeans(bboots, dims = 2)

  irf_mean <- rowMeans(irfboots, dims = 3)
  irf_median <- apply(irfboots, c(1, 2, 3), stats::median)
  irf_lower <- apply(irfboots, c(1, 2, 3), stats::quantile, probs = lower_pctl)
  irf_upper <- apply(irfboots, c(1, 2, 3), stats::quantile, probs = upper_pctl)

  irf_df <- data.frame(
    h = rep(1:n_ahead, each = k * k),
    shock = rep(shock_names, times = k * n_ahead),
    variable = rep(variable_names, each = k, times = n_ahead),
    mean = c(irf_mean),
    median = c(irf_median),
    lower = c(irf_lower),
    upper = c(irf_upper)
  )

  ## Return bootstrap results
  bootstrap <- list(
    VAR = var,
    IRF_df = irf_df,
    IRF_boots = irfboots,
    IRF_mean = irf_mean,
    A_boots = aboots,
    A_mean = a_mean,
    B_boots = bboots,
    B_mean = b_mean,
    id_method = id_method,
    nboot = nboot,
    n_ahead = n_ahead,
    design = design,
    method = method,
    wild_distr = wild_distr,
    bias_adjust = bias_adjust,
    lower_pctl = lower_pctl,
    upper_pctl = upper_pctl
  )

  if (bias_adjust == TRUE) {
    bias <- a_mean - var$A_hat

    var$A_hat <- var$A_hat - bias

    bootstrap <- bootstrap(
      var = var, id_method = id_method, nboot = nboot, n_ahead = n_ahead,
      design = design, method = method, wild_distr = wild_distr,
      bias_adjust = FALSE,
      lower_pctl = lower_pctl, upper_pctl = upper_pctl, ...
    )

    bootstrap$bias_adjust <- TRUE
    bootstrap$bias <- bias
  }

  return(bootstrap)
}
