#' Piece wise polinomial expansion for X (PWP)
#'
#' Calculates and returns a list of matrixes, each one representing a PWP
#' expansion for dimention d. Combines all of the parameters on a relatively
#' fast computation of basis expansion for X described on the thesis and on
#' (Denison, Mallik and Smith, 1998).
#'
#' @inheritParams bpwpm_gibbs
#' @inheritParams calculate_F
#' @param t matrix containing the nodes in which to split the Piecewise Polinomials
#'
#' @return A list of PWP expansion matrixes for each dimention d.
#'
calculate_Phi <- function(X, M, J, K, d, t){

    # Phi is a list containing the basis transformations matrices for each dimension j.
    # For now, this basis expansion is done following the formula on the thesis, optimized as much as posible
    Phi <- list()
    for(j in seq(1,d)){

        # Creating the first regular polinomial
        Phi_partial <- sapply(X = seq(0, M-1), FUN = function(x,y){y^x}, y = X[ , j])

        # Piecewise part
        for(k in seq(1,J-1)){

            # A diagram of this can be found on the thesis
            # Note that in the limit case that K = 0 we need to make adjustments since 0^0 = 1
            if(K != 0){
                Phi_partial <- cbind(Phi_partial, sapply(X = seq(K, M-1),
                                                         FUN = function(x,y){y^x},
                                                         y = pmax(0, X[ ,j] - t[k,j])))
            }
            else{
                temp <- pmax(0, X[ ,j] - t[k,j])
                temp[temp > 0] <- 1
                Phi_partial <- cbind(Phi_partial, temp)

                if(M > 1){
                    Phi_partial <- cbind(Phi_partial, sapply(X = seq(1,M-1),
                                                             FUN = function(x,y){y^x},
                                                             y = pmax(0, X[ ,j] - t[k,j])))
                }
            }
        }

        Phi[[j]] <- Phi_partial
    }

    return(Phi)

}
