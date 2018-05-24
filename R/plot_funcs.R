# Plot Functionality for package bpwpm

plot.bpwpm <- function(object, n = 100, ...){

    if(!('bpwpm' %in% class(object))){
        error("Object not of the class bpwpm")
        geterrmessage()
    }

    # Betas
    p <- plot_chains(object$betas, n, title = "Betas")
    print(p)
    readline(prompt = "Press [enter] to view next plot")
    p <- plot_hist(object$betas, n, title = "Betas")
    print(p)

    for(i in seq(1, length(object$w))){
        readline(prompt = "Press [enter] to view next plot")
        p <- plot_chains(object$w[[i]], n, title = paste("w_",i))
        print(p)
        readline(prompt = "Press [enter] to view next plot")
        p <- plot_hist(object$w[[i]], n, title = paste("w_",i))
        print(p)

    }

}

#-------------------------------------------------------------------------------

#' Plot MCMC Chains
#'
#' Plots the last n draws of an MCMC chain
#'
#' @param mcmc_chain An MCMC Chain matrix. (draws * number of params)
#' @param number The number of draws to plot.
#' @param title Title for the plot
#'
#' @return A ggplot2 plot
#' @export
#'
#' @examples plot_chains(betas), plot_chains(w_j, 1000)
plot_chains <- function(mcmc_chain, number = 100, title = ""){

    dim_mcmc <- dim(mcmc_chain)
    n <- min(number, dim_mcmc[1])

    beta_temp <- tidyr::gather(mcmc_chain[seq(dim_mcmc[1] - n + 1,dim_mcmc[1]),],
                               key = Parameters)

    ggplot2::ggplot(beta_temp, aes(x = rep(seq(1,n),dim_mcmc[2]),
                                   y = value, group = Parameters,
                                   colour = Parameters)) +
             geom_line() + xlab("Index") + ylab("Value") + ggtitle(title)
}

#-------------------------------------------------------------------------------

plot_hist <- function(mcmc_chain, number = 100, title = "", ...){

    dim_mcmc <- dim(mcmc_chain)
    n <- min(number, dim_mcmc[1])

    beta_temp <- tidyr::gather(mcmc_chain[seq(dim_mcmc[1] - n + 1,dim_mcmc[1]),],
                               key = Parameters)

    ggplot2::ggplot(beta_temp, aes(x = value, fill = Parameters)) +
        geom_histogram(..., position = "dodge") + xlab("Value") +
        ggtitle(title)


}

#-------------------------------------------------------------------------------

plot.bpwpm_prediction <- function(object, ...){

    if(!('bpwpm_prediction' %in% class(object))){
        error("Object not of the class bpwpm_prediction")
        geterrmessage()
    }

    plot_each_F(object$Y, object$X, object$bpwpm_params)
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
#' @param F_mat The F matrix calculated via \code{\link{calculate_F}} or
#' alternativly you can pass it the parameters calculated by function
#' \code{\link{posterior_params}}
#'
#' @return d plots for each dimention
#' @export
#'
plot_each_F <- function(Y, X, F_mat){

    if(class(F_mat) == "bpwpm_params"){
        F_mat <- F_mat$estimated_F
    }

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

# Methods for ploting 2D Graphs
#-------------------------------------------------------------------------------


    # Sanitizing Inputs
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

    if(length(Y) != dim(X)[1]){
        error("Y and X have a diferent number of observations")
        geterrmessage()
    }

    # Normal Data
    p <- plot_2D_data(Y,X)
    print(p)
    readline(prompt = "Press [enter] to view next plot")
    p <- plot_2D_proj(Y, X, bpwpm_params, n, alpha)
    print(p)
    readline(prompt = "Press [enter] to view next plot")
    p <- plot_3D_proj(X, bpwpm_params, n, f_of_0)
    print(p)

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

#' Plot 2D projection of the Model
#'
#' Once a model has been run and evaluated, in case that we have a 2D input
#' matrix, we can plot the projection to evaluate the model and its
#' corresponding binary outcomes. Instead of plotting the corresponding conotur
#' of the 3D function ploted by \code{\link{plot_3D_proj}} the output is
#' converted to its corresponding output and mapped to the 2D input space.
#'
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

    if(dim(X)[2] != 2){
        error("Only a 2D plot can be made. X matrix has diferent dimensions")
        geterrmessage()
    }

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
