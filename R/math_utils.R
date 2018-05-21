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

    if(class(Y) == "factor"){
        Y <- as.integer(Y) -1
    }

    p_corrected = pmin(pmax(p, eps), 1-eps)
    ll <- - sum (Y * log(p_corrected) + (1 - Y) * log(1 - p_corrected))/length(Y)

    if(verb){
        cat("\nLog_Loss: ",ll, sep ="")
    }

    return(- sum (Y * log(p_corrected) + (1 - Y) * log(1 - p_corrected))/length(Y))
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

#-------------------------------------------------------------------------------

#' Calculates trained model for new Data
#'
#' Given a set of parameters inherited by the \code{\link{bpwpm_gibbs}} training
#' procedure, we calculate the value of the projection function for new data.
#' This function is both used on the predict and plot_3D functions.
#' @param new_data A new set of data for which to calculate the f(x) projection
#' function
#' @param bpwpm_params A list of bpwpm parameters created by the function
#' \code{\link{posterior_params}}
#'
#' @return A vector of the projection vector for a given set of points
#' @export
#'
#' @examples (test_data, mean_params)
model_projection <- function(new_data, bpwpm_params){

    if(class(bpwpm_params) != "bpwpm_params"){
        error("Invalid bpwpm parameters")
        geterrmessage()
    }

    Phi <- calculate_Phi(X = new_data,
                         M = bpwpm_params$M, J = bpwpm_params$J,
                         K = bpwpm_params$K, d = bpwpm_params$d,
                         t = bpwpm_params$tau)

    F_mat <- calculate_F(Phi = Phi, bpwpm_params$w, d = bpwpm_params$d,
                         intercept = bpwpm_params$intercept)

    return(calculate_projection(F_mat = F_mat, betas = bpwpm_params$betas))
}

#-------------------------------------------------------------------------------

#' Calculate Probabilities of Binary Outcome
#'
#' Given a model, we can calculate the corresponding the probabilites of the
#' random response variable Y. Whereas this is new data or the one used to train
#' the model. Since the model is a probit GLM at this point, we only need to
#' calculate the projection and then plug them on the inverse of the normal
#' gaussian acomulation function
#'
#' @inheritParams model_projection
#'
#' @return a vector of fitted probabilities for the response variable Y
#' @export
#'
posterior_probs <- function(new_data, bpwpm_params){
    if(class(bpwpm_params) != "bpwpm_params"){
        error("Invalid bpwpm parameters")
        geterrmessage()
    }

    z <- model_projection(new_data = new_data, bpwpm_params)

    return(pnorm(z))
}

#-------------------------------------------------------------------------------

#' Calculate Acurracy of the model
#'
#' Given a set of true values and their corresponding fitted probabilities,
#' the function calculates the acurracy of the model defined by:
#'  $1- ,#wrong prediction/# of observations$
#'
#' @inheritParams  log_loss
#'
#' @return The acurracy of the model, given the fitted probabilities and new data
#' @export
#'
#' @examples (new_Y_ata, fitted_probs_for_data)
accurracy <- function(new_Y, p, verb = TRUE){

    if(class(new_Y) == "factor"){
        new_Y <- as.integer(Y) - 1
    }

    n <- length(new_Y)

    predicted_Y <- rep(0, times = n)
    predicted_Y[p > 0.5] <- 1

    wrongs <- sum(abs(new_Y - predicted_Y))

    if(verb) cat(wrongs, " incorrect categorizations\n", sep = "")

    return(1 - wrongs/n)
}


