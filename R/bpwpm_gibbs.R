#' Bayesian Piecewise Polinomial Model
#'
#' Gibbs Sampler for simulating draws of parameters for the Bayesian Piece Wise
#' Polinomial model described on the thesis.
#'
#' @param Y Response vector of n binary observatios (integers 0,1 - vector of
#'   size n)
#' @param X Design matrix of n observations and d covariables (numeric - n*d)
#' @param M M minus 1 is the degree of the polinomial (integer - M > 0)
#' @param J Númber of intervals in each dimention (integer - J > 1)
#' @param K Order of continuity in the derivatives (integrer - K < M)
#' @param intercept Intercept:= To add an itercept term to the estimation. It's
#'   recomended to have one (logical)
#' @param precision_beta If using the default sigmas for betas, a diagonal
#'   matrix will be used. Precision controls its magnitude (numeric - precision
#'   > 0)
#' @param precision_w If using the default sigmas for w, a diagonal matrix will
#'   be used. Precision controls its magnitude (numeric - precision > 0)
#' @param draws Númber of samples to draw from the Gibbs Sampler (integer - draw
#'   > 0)
#' @param z_init z_init:= Initial values for Auxiliary variable z (numeric -
#'   vector of size n)
#' @param beta_init Initial value for the Gibbs Sampler Chain (numeric - vector
#'   of size d)
#' @param mu_beta_0 Prior Mean of Beta (numeric - vector of size d)
#' @param sigma_beta_0_inv Prior Inverse of the Variance-Covariance Matrix of
#'   beta (numeric - matrix of size d*d)
#' @param w_init Inital value for the Gibbs Sampler Chain (numeric - matrix of
#'   size N*d)
#' @param mu_w_0 Prior Mean of Beta (numeric - matrix of size N*d)
#' @param sigma_w_0_inv sigma_w_0_inv:= Prior Inverse if the Variance-Covariance
#'   Matrices of w (list - d elements, each element is a numeric matrix of size
#'   N*N)
#' @param verb verbose, if True, prints aditional information (logical)
#' @param debug If TRUE, print even more info to help with debugging (logical)
#' @param keep If TRUE, saves a copy of the transformations Phi (logical)
#'
#' @return A list containing all the elements of the model
#' @export
#'
#' @examples See the main document of the thesis for a couple of full examples
#' with its corresponding analysis.

bpwpm_gibbs <- function(Y, X, M, J, K,
                  intercept = TRUE, precision_beta = 1, precision_w = 1,
                  draws = 10^3,
                  z_init = NULL, beta_init = NULL, mu_beta_0 = NULL,
                    sigma_beta_0_inv = NULL,
                  w_init = NULL, mu_w_0 = NULL, sigma_w_0_inv = NULL,
                  verb = TRUE, debug = FALSE, keep = TRUE){

    # 0. Stage Setting----------------------------------------------------------

    # 0.1 initializing Data
    dim_X <- dim(X)
    n <- dim_X[1] # Number of observations
    d <- dim_X[2] # Number of covariables

    N <- M * J - K * (J - 1) # Númber of basis funcitons

    # Z Initial parameters
    if(is.null(z_init)){
        z <- rep(0,n) # Standard z
    }else{
        z <- z_init
    }

    # Initial parameters for beta
    if(intercept){
        if(is.null(beta_init)){
            betas <- rep(0,d+1) # Standard beta
        }else{
            betas <- beta_init
        }
        if(is.null(mu_beta_0)){
            mu_beta_0 <- rep(0,d+1)
        }
        if(is.null(sigma_beta_0_inv)){
            sigma_beta_0_inv <- precision_beta*diag(d+1)
        }

    }else{
        if(is.null(beta_init)){
            betas <- rep(0,d) # Standard beta
        }else{
            betas <- beta_init
        }
        if(is.null(mu_beta_0)){
            mu_beta_0 <- rep(0,d)
        }
        if(is.null(sigma_beta_0_inv)){
            sigma_beta_0_inv <- precision_beta*diag(d)
        }
    }

    # Initial parameters for w
    if(is.null(w_init)){
        # Standard w, the weight are taken uniformly throught w
        w <- matrix(1/N, nrow = N, ncol = d)
    }else{
        w <- w_init
    }

    if(is.null(mu_w_0)){
        mu_w_0 <- matrix(0L, nrow = N, ncol = d)
    }

    if(is.null(sigma_w_0_inv)){
        sigma_w_0_inv <- rep(list(precision_w*diag(N)),d)
    }

    # O.2 Dimension check
    if(!all(identical(dim_X[1],length(Y)), K < M, M > 0, J > 1)){
        stop("Error in dimensionalities")
        geterrmessage()
    }

    rm(dim_X)

    # 0.3 Basic Info print
    cat("\tBPWPM MODEL\n\t
        Dimensions and pararamters check\n\t",
        "Algorithm based on d = ", d, " covariables", "\n\t",
        "Number of nodes J - 1 = ", J - 1, "\n\t",
        "Order of polinomial M - 1 = ", M - 1, "\n\t",
        "Order of continuity on derivatives K = ", K, "\n\t",
        "Number of basis functions N = ", N, "\n\t",
        "Number of observations n = ", n, "\n\n", sep = "")

    # 1. Initializing Gibbs Sampler---------------------------------------------

    cat("Initializing Gibbs Sampler\n\n")

    # 1.1 Node Initialization.
    # Setting Nodes on the quantiles of X. (numeric, (J-1)*d)
    t <- sapply(X, quantile, probs = seq(0,1, by = 1/J))
    t <- matrix(t[-c(1,J+1), ], nrow = J-1, ncol = d)

    if(verb){
        # Nodes
        cat("\tNodes Locations\n")
        print(t, digits = 3)
        cat("\n")

        # Weights
        cat("\nInitial Weights\n")
        print(w, digits = 3)
        cat("\n")
    }

    # 1.2 Piece Wise Polinomial expansion Phi
    Phi <- calculate_Phi(X = X, M = M, J = J, K = K, d = d, t = t)

    # 1.3 Calculating initial F
    F_mat <- calculate_F(Phi = Phi, w = w, d = d, intercept = intercep)

    if(verb){
        cat("\tInitial F\n")
        print(F_mat, digits = 3)
        cat("\n")
    }

    # 1.4 Final variables
    # List for storing w.
    sim_w <- list()

    # For speed calculation
    Phi_cross <- lapply(X = Phi, FUN = function(x){crossprod(x,x)})

    # 2. Gibbs Sampler, z -> beta -> w -> F -> z -> ... ------------------------

    for(k in seq(1,draws)){

        cat("\nIter: ", k, sep ="")

        # Linear predictor (Auxiliary variable)
        eta <- crossprod(t(F_mat),betas)

        # 2.1. Z - Sampling from the truncated normal distribution for Z.
        z[Y == 0] <- qnorm(runif(n = sum(1-Y),0,
                                 pnorm(0,eta[Y == 0],1)),eta[Y==0],1)
        z[Y == 1] <- qnorm(runif(sum(Y),
                                 pnorm(0,eta[Y == 1],1),1), eta[Y == 1], 1)


        # 2.2. BETA - Sampling from the final distribution for beta
        sigma_beta <- solve(sigma_beta_0_inv + crossprod(F_mat,F_mat))
        mu_beta <- crossprod(t(sigma_beta),(crossprod(t(sigma_beta_0_inv),mu_beta_0) + crossprod(F_mat,z)))

        betas <- c(rmvnorm(1,mu_beta,sigma_beta)) # Simulating from the resulting distribution

        # Making beta simulation matrix
        if(k == 1){
            sim_beta <- betas
        }else{
            sim_beta <- rbind(sim_beta,betas)
        }

        if(verb){
            cat("\n\tBeta:\t", format(betas,digits = 2, width = 10), sep = "")
        }

        # 2.3. W - Sampling from the distribution for w
        for(j in 1:d){

            # Residuals for dim j
            h <- (z - crossprod(t(F_mat[, -( j+1 )]), betas[-(j+1)]))/betas[j+1]

            if(debug){cat("\n\tPartial Residuals:", format(head(h),digits = 2, width = 10), cat = "")}

            sigma_w_j <- solve(sigma_w_0_inv[[j]] + Phi_cross[[j]])
            mu_w_j <- crossprod(t(sigma_w_j),(crossprod(t(sigma_w_0_inv[[j]]),mu_w_0[ ,j]) + crossprod(Phi[[j]],h)))

            if(debug){cat("\n\tmu_w", j,":\t", format(mu_w_j, digits = 2, width = 10), cat = "")}

            w_j <- c(rmvnorm(1,mu_w_j, sigma_w_j))

            if(verb){cat("\n\tw",j,":\t", format(w_j, digits = 2, width = 10),sep = "")}
            if(debug){cat("\n\tF", j,"previous update:\t", format(head(F_mat[,(j+1)]), digits = 2, width = 10), cat = "")}

            F_mat[ , (j+1)] <- crossprod(t(Phi[[j]]), w_j)

            if(debug){cat("\n\tF", j,"after update:\t", format(head(F_mat[,(j+1)]), digits = 2, width = 10), cat = "")}
            if(debug){cat("\n\tF", 1,":\t", format(head(F_mat[,1]), digits = 2, width = 10), cat = "")}

            # Making sim_w
            if(k == 1){
                sim_w[[j]] <- w_j
            }else{
                sim_w[[j]] <- rbind(sim_w[[j]],w_j)
            }
        }
    }

    # 3. Creating Output--------------------------------------------------------

    #

    if(keep){
        model <- list(betas = sim_beta, w = sim_w, F_mat = F_mat, Phi = Phi)
    } else{
        model <- list(betas = sim_beta, w = sim_w, F_mat = F_mat)
    }

    return(model)
}
