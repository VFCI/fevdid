#' Convert a VAR with p lags to a representation with just one lag
#'
#' @param svar svar object, from id.chol
#'
#' @return matrix coefficients for VAR 1 representation
#'
svar_to_svar1 <- function(svar) {

    ## Check that the VAR is an svar
    if (!inherits(svar, "svars")) stop("Please pass a svars.")

    ## Get coeefficient matrix, ignore constant
    if (svar$type == "const") {
        a <- svar$A_hat[, -1]
    } else {
        a <- svar$A_hat
    }

    ## Number of variables
    k <- dim(svar$y)[2]

    ## Number of lags
    p <- svar$p

    ## Create Var(1) objects
    my <- cbind(diag(k * (p - 1)), matrix(0, k * (p - 1), k))
    mx <- rbind(a, my)
    me <- rbind(svar$B, matrix(0, k * (p - 1), k))

    return(list(my = my, mx = mx, me = me))
}
