#' Bayesian Piecewise Polinomial Model
#'
#' Gibbs Sampler for simulating draws of parameters for the Bayesian Piece Wise
#' Polinomial model described on the thesis.
#'
#' @param Y Response vector of n binary observatios (integers 0,1 - vector of
#'   size n) Can be encoded as a factor a numeric vector.
#' @param X Design matrix of n observations and d covariables (numeric - n*d)
#' @param M M minus 1 is the degree of the polinomial (integer - M > 0)
#' @param J Number of intervals in each dimention (integer - J > 1)
#' @param K Order of continuity in the derivatives (integrer - K < M)
#' @param precision_w If using the default sigmas for w, a diagonal matrix will
#'   be used. Precision controls its magnitude (numeric - precision > 0)
#' @param draws Númber of samples to draw from the Gibbs Sampler (integer - draw
#'   > 0)
#' @param tau the initial position of nodes selected by the user. although·
#' arbitraty they need to match the dimentions. (numeric - (J-1)*d)
# @param beta_init Initial value for the Gibbs Sampler Chain (numeric - vector
#'   of size d)
# @param mu_beta_0 Prior Mean of Beta (numeric - vector of size d)
# @param sigma_beta_0_inv Prior Inverse of the Variance-Covariance Matrix of
#'   beta (numeric - matrix of size d*d)
#' @param w_init Inital value for the Gibbs Sampler Chain (numeric - matrix of
#'   size N*d)
#' @param mu_w_0 Prior Mean of Beta (numeric - matrix of size N*d)
#' @param sigma_w_0_inv sigma_w_0_inv:= Prior Inverse if the Variance-Covariance
#'   Matrices of w (list - d elements, each element is a numeric matrix of size
#'   N*N)
#' @param indep_terms Keeps the independent terms in the PWP expansion.
#'   Leaving those terms might lead to identificability errors. (logical)
#' @param verb short for verbose, if TRUE, prints aditional information (logical)
#' @param debug If TRUE, print even more info to help with debugging (logical)
#'
#' @return An object of the class "bpwpm" containing at the following
#' components:
#' \describe{
#' \item{betas: }{A data frame containing the Gibbs sampler simulation for beta}
#' \item{w: }{A list of d elements. Each one is a data frame containign the
#' simulation of the w_j parameters for each dimetnion j.}
#' \item{Phi: }{The PWP Expansion for input matrix X and nodes selected on
#' percentiles}
#' \item{tau: }{Nodes used for training}
#' \item{M: }{Initial parameters}
#' \item{J: }{Initial parameters}
#' \item{K: }{Initial parameters}
#' \item{d: }{Number of dimentions}
#' \item{indep_terms: }{Logical. If independent terms are keept}
#' \item{info}{A string that prints the basic information of the mode. Used for
#' the summary function.}
#' }
#' @export
#'
#' @examples See the main document of the thesis for a couple of full examples
#' with its corresponding analysis.

bpwpm_gibbs <- function(Y, X, M, J, K,
                draws = 10^3,
                indep_terms = FALSE,
                tau = NULL,
                w_init = NULL,
                mu_w_0 = NULL, sigma_w_0_inv = NULL, precision_w = 1,
                verb = FALSE, debug = FALSE){

    # 0. Stage Setting----------------------------------------------------------

    # 0.1 initializing Data
    dim_X <- dim(X)
    n <- dim_X[1] # Number of observations
    d <- dim_X[2] # Number of covariables

    N <- M * J - K * (J - 1) # Númber of basis funcitons

     # Initial parameters for betas
    betas <- rep(1,d+1)

    # Initial parameters for w
    if(is.null(w_init)){
        # Standard w, the weight are taken uniformly throught w
        if(indep_terms == FALSE){
            N <- N - 1
            w <- matrix(1/N, nrow = N, ncol = d) # Removes independent parameters
            cat("Removing independent terms, using N-1 parameters")
        }else{
            w <- matrix(1/N, nrow = N, ncol = d) # Removes independent parameters
            cat("Independent terms kept, using N parameters")
        }

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

    # Encoding Y
    if(class(Y) == "factor"){
        Y <- as.integer(Y) - 1
    }

    # 0.3 Basic Info print
    info <- paste("\n\tBPWPM MODEL\n\t
    \tDimensions and pararamters check\n\t",
    "Algorithm based on d = ", d, " covariables", "\n\t",
    "Number of nodes J - 1 = ", J - 1, "\n\t",
    "Order of polinomial M - 1 = ", M - 1, "\n\t",
    "Order of continuity on derivatives K = ", K, "\n\t",
    "Number of basis functions N = ", N, "\n\t",
    "Number of observations n = ", n, "\n\n", sep = "")
    cat(info)

    # 1. Initializing Gibbs Sampler---------------------------------------------

    cat("Initializing Gibbs Sampler\n")

    # 1.1 Node Initialization.
    # Setting Nodes on the quantiles of X. (numeric, (J-1)*d)
    if(is.null(tau)){
        tau <- sapply(data.frame(X), quantile, probs = seq(0,1, by = 1/J))
        tau <- matrix(tau[-c(1,J+1), ], nrow = J - 1, ncol = d)
    }else if(!(dim(tau)[1] == (J-1) && dim(tau)[2] == d)){
        error("Dimentions of the given tau does not match up. The standard tau matrix is recomended")
        geterrmessage()
    }

    if(verb){
        # Nodes
        cat("\tNodes Locations\n")
        print(tau, digits = 3)
        cat("\n")

        # Weights
        cat("\tInitial Weights\n")
        print(w, digits = 3)
        cat("\n")
    }

    # 1.2 Piece Wise Polinomial expansion Phi
    Phi <- calculate_Phi(X, M, J, K, d, tau, indep_terms)

    # 1.3 Calculating initial F
    F_mat <- calculate_F(Phi, w, d)

    if(verb){
        cat("\tInitial F\n")
        print(F_mat, digits = 3)
        cat("\n")
    }

    # 1.4 Final variables
    # List for storing w.
    sim_w <- list()

    # For speedy calculation
    Phi_cross <- lapply(X = Phi, FUN = function(x){crossprod(x,x)})
    z <- rep(0, times = n)

    # 2. Gibbs Sampler, z -> beta -> w -> F -> z -> ... ------------------------

    # cat("New Version!")
    for(k in seq(1,draws)){

        # Print number of iterations
        if(verb){cat("\nIter: ", k, sep ="")
        }else if((k %% 100) == 0){
            cat("\nIter: ", k, sep = "")
        }

        # Linear predictor (Auxiliary variable)
        eta <- crossprod(t(F_mat),betas)

        # 2.1. Z - Sampling from the truncated normal distribution for Z.
        eps <- 1e-15

        temp_probs <- pnorm(0,eta,1)
        temp_probs = pmin(pmax(temp_probs, eps), 1-eps)

        z[Y == 0] <- qnorm(runif(n = sum(1-Y), min = 0,
                                 max = temp_probs[Y == 0]), eta[Y == 0],1)
        z[Y == 1] <- qnorm(runif(n = sum(Y),
                                 min = temp_probs[Y == 1],max = 1), eta[Y == 1], 1)

        # 2.2. BETA - Sampling from the final distribution for beta
        #sigma_beta <- solve(sigma_beta_0_inv + crossprod(F_mat,F_mat))
        #mu_beta <- crossprod(t(sigma_beta),(crossprod(t(sigma_beta_0_inv),mu_beta_0) +
        #                                        crossprod(F_mat,z)))
        # betas <- c(mvtnorm::rmvnorm(1,mu_beta,sigma_beta)) # Simulating from the resulting distribution

        # Making beta simulation matrix
        if(k == 1){
            sim_beta <- betas # First iteration
        }else{
            sim_beta <- rbind(sim_beta,betas) # Apending Betas
        }

        # GAM first step
        betas[1] <- mean(eta)

        if(verb){
            cat("\n\tBeta:\t", format(betas,digits = 2, width = 10), sep = "")
        }


        # 2.3. W - Sampling from the distribution for w
        for(j in 1:d){

            # Residuals for dim j
            h <- (eta - crossprod(t(F_mat[, -( j+1 )]), betas[-(j+1)]))/betas[j+1]

            if(debug){cat("\n\tPartial Residuals:", format(head(h),digits = 2, width = 10), cat = "")}

            sigma_w_j <- solve(sigma_w_0_inv[[j]] + Phi_cross[[j]])
            mu_w_j <- crossprod(t(sigma_w_j),(crossprod(t(sigma_w_0_inv[[j]]),mu_w_0[ ,j])
                                              + crossprod(Phi[[j]],h)))

            if(debug){cat("\n\tmu_w", j,":\t", format(mu_w_j, digits = 2, width = 10), cat = "")}

            w_j <- c(mvtnorm::rmvnorm(1,mu_w_j, sigma_w_j))

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

    # 3.1 Naming Betas
    rownames(sim_beta) <- NULL
    colnames(sim_beta) <- paste("beta_", seq(0,d), sep = "")

    # 3.2 Naming and reformating W's
    names(sim_w) <- paste("w_", seq(1,d), sep = "")
    lapply(seq(1,length(sim_w)), function(x){sim_w[[x]] <<- data.frame(sim_w[[x]])})
    lapply(seq(1,length(sim_w)), function(x){row.names(sim_w[[x]]) <<- NULL})
    lapply(seq(1,length(sim_w)), function(x){colnames(sim_w[[x]]) <<-
        paste("w_",x,"_",seq(1,N), sep = "")})

    model <- list(betas = data.frame(sim_beta),
                  w = sim_w,
                  Phi = Phi,
                  tau = tau,
                  M = M,
                  J = J,
                  K = K,
                  d = d,
                  indep_terms = indep_terms,
                  info = info)
    class(model) <- 'bpwpm'

    return(model)
}
