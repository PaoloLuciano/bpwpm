# Math Utils

#-------------------------------------------------------------------------------

#' Log Loss
#'
#' An implementation of the Log-Loss function for the binomial case
#'
#' @inheritParams bpwpm_gibbs
#' @param p Vector of fitted probabilities for each case
#' @param eps Machine error to hanlde limit cases on the logarithm function
#'
#' @return The value of the Log Loss function (numeric). The smaller, the better
#'
#' @examples log_loss(true_values, fitted probabilities)
#' log_loss(true_values, fitted probabilities)
#' #' log_loss(true_values, fitted probabilities, 1e-30, FALSE)
log_loss <- function(Y, p, eps = 1e-15, verb = TRUE){

    p_corrected = pmin(pmax(p, eps), 1-eps)
    ll <- - sum (y * log(p_corrected) + (1 - y) * log(1 - p_corrected))/length(y)

    if(verb){
        cat("\nLog_Loss: ",ll, sep ="")
    }

    return(- sum (y * log(p_corrected) + (1 - y) * log(1 - p_corrected))/length(y))
}

#-------------------------------------------------------------------------------

#' Mode
#'
#' Calculates the mode of a vector. Ties are resolved by the first element
#'
#' @param x A numeric vector
#'
#' @return The mode of the vector
#'
mode <- function(x) {

    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

#-------------------------------------------------------------------------------

#' Calculate Projection Vector
#'
#' @param F_mat The PWP transformed input space, calculated by a set of \code{w}
#'  and a \code{Phi} matrix by the function \code{\link{calculate_F}}
#' @param betas The posterior puntual estimation for beta parameters calculated
#' by \code{\link{posterior_params}}
#'
#' @return A numeric vector representing the projection of \code{R^d} into
#' \code{R} given all of the parametres
#' @export
#'
calculate_projection <- function(F_mat, betas){
    return(crossprod(t(F_mat),betas))
}
