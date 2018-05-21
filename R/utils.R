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

    params <- list(betas = estimated_betas, w = estimated_w,
                   tau = bpwpm_model$tau,
                   M = bpwpm_model$M,
                   J = bpwpm_model$J,
                   K = bpwpm_model$K,
                   d = bpwpm_model$d,
                   intercept = bpwpm_model$intercept)
    class(params) <- 'bpwpm_params'

    return(params)
}

#-------------------------------------------------------------------------------

#' Plot MCMC Chains
#'
#' Plots the last n draws of an MCMC chain
#'
#' @param mcmc_chain An MCMC Chain matrix. (draws * number of params)
#' @param number The number of draws to plot.
#'
#' @return A ggplot2 plot
#' @export
#'
#' @examples plot_chains(betas), plot_chains(w_j, 1000)
plot_chains <- function(mcmc_chain, number = 100){

    dim_mcmc <- dim(mcmc_chain)
    n <- min(number, dim_mcmc[1])

    beta_temp <- tidyr::gather(mcmc_chain[seq(dim_mcmc[1] - n + 1,dim_mcmc[1]),],
                               key = Parameters)

    ggplot2::ggplot(beta_temp, aes(x = rep(seq(1,n),dim_mcmc[2]),
                                   y = value, group = Parameters, colour = Parameters)) +
        geom_line() + xlab("Index") + ylab("Value")
}

#-------------------------------------------------------------------------------

#' Plots each F(X)
#'
#' With the \code{w} parameters calculated from the Ginns run, and
#' \code{\link{posterior_params}}, a Final F matrix can be calculated. and
#' hence, ploted against every Input X to see how does the PWP expansion looks
#' like for the specified set of parameters.
#'
#' @param Y A vector of binary response. Can be encoded as either a factor
#' vector or as a numeric one.
#' @param X A data frame or matrix containing the original Inputs for the model.
#' @param F_mat The F matrix calculated via \code{\link{calculate.F}}
#'
#' @return d plots for each dimention
#' @export
#'
plot_each_F <- function(Y,X,F_mat){

    d <- dim(X)[2]

    if(class(Y) == "numeric" | class(Y) == "integer"){
        Y <- as.factor(Y)
    }

    # If F has an intercept term
    if(dim(F_mat)[2] == (d + 1)){
        F_mat <- F_mat[, -1]
    }

    for(i in seq(1:d)){
        p <- ggplot2::qplot(x = X[,i],  y = F_mat[,i],
                                    color = Y) +
                                    xlab(paste("X_",i, sep = "")) +
                                    ylab(paste("F_",i,"(X_",i,")",sep = ""))
        print(p)
        if(i != d){
            readline(prompt="Press [enter] to view next plot")
        }
    }
}

#-------------------------------------------------------------------------------

#' Scatter Plot of 2D data
#'
#' @inheritParams plot_each_F
#' @param X Input Matrix of 2D (2 columns).
#'
#' @return A scatter plot of the 2 groups
#' @export
#'
#' @examples (Y = rbinom(100, 1, 4), X = cbind(rnorm(100), rnorm(100)))
plot_2D_data <- function(Y,X){

    if(dim(X)[2] != 2){
        error("Only a 2D plot can be made. X matrix has diferent dimensions")
        geterrmessage()
    }

    if(class(Y) == "numeric"){
        Y <- factor(Y)
    }

    if(class(X) == "matrix"){
        X <- data.frame(X)
    }

    ggplot2::ggplot(data = X, aes(x = X[, 1], y = X[ ,2], col = Y)) +
        geom_point() + xlab("X_1") + ylab("X_2")
}

#-------------------------------------------------------------------------------

#' Plots the 3D representation of the projection function
#'
#' Given the set of parmeters and the input data in 2D, this function calculates
#' and plots the data on a 3D linear space defined by the input matrix X.
#' @inheritParams plot_2D_data
#' @inheritParams model_projection
#' @param n How thin is the grid to be made
#' @param f_of_0 If the constant function 0 is to be ploted
#'
#' @return a 3d WireFrame Plot
#' @export
#'
plot_3D_proj <- function(X, bpwpm_params, n, f_of_0 = TRUE){

    mins <- apply(X, 2, min)
    maxs <- apply(X, 2, max)

    linspace <- expand.grid(X1 = seq(mins[1], maxs[1], by = 1/n),
                            X2 = seq(mins[2], maxs[2], by = 1/n))

    linspace$f <- model_projection(new_data = linspace,
                                   bpwpm_params = bpwpm_params)

    if(f_of_0){

        m <- dim(linspace)[1]
        linspace$g <- rep(1, times = m)

        # Building the 0 gridspace
        linspace <- rbind(linspace, data.frame(X1 = linspace[,1], X2 = linspace[,2]
                                               ,f = rep(0, times = m),
                                               g = rep(0, times = m)))
    }

    lattice::wireframe(f ~ X1 * X2, data = linspace, group = g,
                       drape = TRUE,
                       aspect = c(1,1),
                       main = paste("3D plot for: M = ", bpwpm_params$M,
                                    ", J = ", bpwpm_params$J, ", K = ", bpwpm_params$K),
                       frame.plot = FALSE,
                       colorkey = FALSE,
                       scales = list(arrows = FALSE))
                       # col.groups = rgb(c(255,0,0), c(0,255,0), alpha = 70,maxColorValue = 255),
                       # col.regions = colorRampPalette(c("blue", "red"))(50))
                       # at = 0, col.regions = c("red", "blue"))

}

#-------------------------------------------------------------------------------

#' Plot 2D projection of the Model
#'
#' Once a model has been run and evaluated, in case that we have a 2D input
#' matrix, we can plot the projection to evaluate the model and its
#' corresponding binary outcomes. Instead of plotting the corresponding conotur
#' of the 3D function ploted by \link{\code{plot_3D_proj}} the output is
#' converted to its corresponding output and mapped to the 2D input space.
#'
#' @inheritParams plot_2D_data
#' @inheritParams plot_3D_proj
#' @param alpha the corresponding alpha transparency param for the output space.
#'
#' @return A ggplot2 scatter plot
#' @export
#'
plot_2D_proj <- function(Y, X, bpwpm_params, n, alpha = 0.6){

    if(dim(X)[2] != 2){
        error("Only a 2D plot can be made. X matrix has diferent dimensions")
        geterrmessage()
    }

    if(class(Y) == "numeric"){
        Y <- factor(Y)
    }

    if(class(X) == "matrix"){
        X <- data.frame(X)
    }

    mins <- apply(X, 2, min)
    maxs <- apply(X, 2, max)

    linspace <- expand.grid(X1 = seq(mins[1] - 0.2, maxs[1] + 0.2, by = 1/n),
                            X2 = seq(mins[2] - 0.2, maxs[2] + 0.2, by = 1/n))

    linspace$Y <- model_projection(new_data = linspace,
                                   bpwpm_params = bpwpm_params)

    linspace$Y<-  as.factor(as.integer(linspace$Y >= 0))
    linspace$a <- rep(alpha, times = dim(linspace)[1])

    data <- data.frame(cbind(X,Y), a = rep(1, times = dim(X)[1]))
    colnames(data) <- c("X1","X2","Y","a")

    data <- data.frame(rbind(data,linspace))

    ggplot2::ggplot(data = data) +
    geom_point(aes(x = X1, y = X2, col = Y, alpha = a), show.legend = FALSE) +
        xlab("X_1") + ylab("X_2")

}

#-------------------------------------------------------------------------------


