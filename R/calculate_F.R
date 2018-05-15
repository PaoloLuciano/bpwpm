#' F matrix calculation
#'
#' Private function for calculating the F matrix, described on the thesis.
#' This is the transformed input matrix that depends on the piecewise polinomial
#' expansion Phi and a set of weights w.
#'
#' @param Phi Piecewise Polinomail expansion for an input matrix X previosly
#'  calculated by \code{\link{calculate_Phi.R}}
#' @param w Set of weights for which to calculate F
#' @param d Number of dimentions, this parameter helps to improve efficiency
#' @param intercept Intercept optional paramter. Logical
#'
#' @return F matrix
#'
calculate_F <- function(Phi, w, d, intercept){


    # Calculating F matrix
    mat_F <- crossprod(t(Phi[[1]]),w[,1])

    if(d>1){
        for(j in seq(2,d)){
            mat_F <- cbind(mat_F,crossprod(t(Phi[[j]]),w[,j]))
        }
    }

    # Adding intercept terms
    if(intercept){
        mat_F <- cbind(rep(1,dim(mat_F)[1]), mat_F)
    }
    return(mat_F)
}
