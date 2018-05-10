calculate_F <- function(Phi, w, d, intercept = TRUE){


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
