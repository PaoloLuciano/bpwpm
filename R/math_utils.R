# Math Utils

#' Log Loss
#'
#' An implementation of the Log-Loss function for the binomial case
#'
#' @inheritParams bpwpm_gibbs
#' @param p Vector of fitted probabilities for each case
#' @param eps Machine error to hanlde limit cases on the logarithm function
#'
#' @return The value of the Log Loss function
#'
#' @examples log_loss(true_values, fitted probabilities)
#' log_loss(true_values, fitted probabilities)
#' #' log_loss(true_values, fitted probabilities, 1e-30, FALSE)
log_loss <- function(Y, p, eps = 1e-15, verb = TRUE){

    "Function that evaluates if the log loss of a binary probability model is withina a certain threshold

    INPUTS:
    y:= Vector of observed binary outcomes (integer - n)
    p:= Vector of fitted probabilities ie: p = P(y = 1) (numeric - n)
    eps:= Error for evaluating the logarithm function. Important for numerical stability (numeric)
    verb:= Boolean for printing information (boolen)

    OUTPUT:
    LogLoss:= value of the LogLoss function for binary outcomes (numeric)
    "
    p_corrected = pmin(pmax(p, eps), 1-eps)
    ll <- - sum (y * log(p_corrected) + (1 - y) * log(1 - p_corrected))/length(y)

    if(verb){
        cat("\nLog_Loss: ",ll, sep ="")
    }

    return(- sum (y * log(p_corrected) + (1 - y) * log(1 - p_corrected))/length(y))
}
