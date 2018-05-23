# Series of utility functions to evaluate the MCMC chain for the bpwpm model

#-------------------------------------------------------------------------------

#' Thin MCMC Chain
#'
#' Thins and eliminates the burn in period of the Gibbs Chain.
#'
#' @param mcmc_chain An MCMC Chain matrix. (draws * number of params)
#' @param thin A thinning parameter for the MCMC Chain
#' @param burn_in A burn in parameter for the MCMC Chain
#'
#' @return The thinned down version of the MCMC chain
#' @export
#'
#' @examples (beta, 2, 2000)
thin_chain <- function(mcmc_chain, thin = 2, burn_in = 2000){

    draws <- dim(mcmc_chain)[1]

    if(burn_in > draws){stop("Burn In parameter too large")
        geterrmessage()}

    return(mcmc_chain[seq(1, floor((draws - burn_in)/(thin + 1))) + burn_in, ])
}

#-------------------------------------------------------------------------------

#' Thinning of a BPWPM
#'
#' Thins all of the MCMC chains of a bpwpm and returns an object of the same kind
#'
#' @param bpwpm_model A bpwpm object created by \code{\link{bpwpm_gibbs}}
#' @inheritParams thin_chain
#'
#' @return An object of the same kind of the bpwpm but thinned down
#' @export
#'
thin_bpwpm <- function(bpwpm_model, thin = 2, burn_in = 2000){

    if(class(bpwpm_model) != 'bpwpm'){
        error("Object not of class bpwpm")
        geterrmessage()
    }

    bpwpm_model_copy <- bpwpm_model
    bpwpm_model_copy$betas <- thin_chain(bpwpm_model$betas, thin = thin, burn_in = burn_in)
    bpwpm_model_copy$w <- lapply(bpwpm_model$w, thin_chain, thin = thin, burn_in = burn_in)

    return(bpwpm_model_copy)
}

#-------------------------------------------------------------------------------

#' Puntual estimation for parameters
#'
#' Given a model output by the function bpwpm_gibbs, it thins down the chain and
#' makes puntual estimation for the parameters. This parameter structure can
#' later be used to calculate the F matrix and other results.
#'
#' @inheritParams thin_bpwpm
#' @param type The type of punctual estimation for the parameters. Options
#' include: mean, mode or median.
#'
#' @return A parameter list that contains the estimations of beta and w.
#' @export
#'
#' @examples (model, 2, 2000, 'mean') (model, 0, 0, 'median')
posterior_params <- function(bpwpm_model, thin, burn_in, type = 'mean'){

    if(class(bpwpm_model) != 'bpwpm'){
        error("Object not of class bpwpm")
        geterrmessage()
    }

    if(type == 'mean'){func <- mean}
    else if(type == 'mode'){func <- mode}
    else if(type == 'media'){func <- median}
    else{error("Incorrect type of parameter estimation")
        geterrmessage()}

    # Fist the model is thinned down
    thined_model <- thin_bpwpm(bpwpm_model, thin = thin, burn_in = burn_in)

    # Estimation for beta
    estimated_betas <- sapply(thined_model$betas, func)

    # Estimation for w
    estimated_w <- lapply(seq(1:length(thined_model$w)), function(x){
        sapply(thined_model$w[[x]], func)})
    estimated_w <- data.frame(matrix(unlist(estimated_w), ncol = length(estimated_w)))
    colnames(estimated_w) <- paste("w_", seq(1,dim(estimated_w)[2]), sep = "")

    estimated_F <- calculate_F(bpwpm_model$Phi,estimated_w, d = bpwpm_model$d,
                               intercept = bpwpm_model$intercept)

    params <- list(betas = estimated_betas, w = estimated_w,
                   tau = bpwpm_model$tau,
                   estimated_F = estimated_F,
                   M = bpwpm_model$M,
                   J = bpwpm_model$J,
                   K = bpwpm_model$K,
                   d = bpwpm_model$d,
                   intercept = bpwpm_model$intercept
                   )
    class(params) <- 'bpwpm_params'

    return(params)
}
