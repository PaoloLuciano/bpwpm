summary.bpwpm <- function(object,
                          digits = max(3L, getOption("digits") - 3L), ...) {

    if(!('bpwpm' %in% class(object))){
        error("Object not of the class bpwpm")
        geterrmessage()
    }

    cat(object$info)

    print_info <- function(object, title, digits){
        cat("\n",title,"\n")
        qq <- sapply(object, quantile)
        temp <- do.call(data.frame,
                        list(min = signif(qq[1, ], digits),
                             FirstQ. = signif(qq[2, ], digits),
                             mean = signif(sapply(object, mean), digits),
                             Median = signif(qq[3, ], digits),
                             ThirdQ. = signif(qq[4, ], digits),
                             Max = signif(qq[5, ], digits),
                             Sd. = signif(sapply(object, sd), digits),
                             Mode = signif(sapply(object, mode), digits)))
        print(t(temp))
    }

    print_info(object$betas, "Betas", digits)

    for(i in seq(1, length(object$w))){
        print_info(object$w[[i]], paste("w_",i, sep = ""), digits)
    }
}

#-------------------------------------------------------------------------------

summary.bpwpm_params <- function(object,
                                 digits = max(3L, getOption("digits") - 3L),
                                 verb = FALSE, ...){

    if(!('bpwpm_params' %in% class(object))){
        error("Object not of the class bpwpm_params")
        geterrmessage()
    }

    cat("\nPosterior Estimated Params: \n")

    print(object$betas)
    cat("\n")
    print(object$w)

    if(verb){
        cat("\nNodes\n")
        print(object$tau)

        # cat("\nF_Transformation\n")
        # print(object$params$estimated_F)
    }
}

#-------------------------------------------------------------------------------

summary.bpwpm_prediction <- function(object,
                                    digits = max(3L, getOption("digits") - 3L),
                                    verb = FALSE, ...){

    if(!('bpwpm_prediction' %in% class(object))){
        error("Object not of the class bpwpm_prediction")
        geterrmessage()
    }

    cat(object$info)
    cat("Prediction Results:\n",
        "\tAccuracy:\t", signif(object$accuracy*100, (digits - 1)), "%\n",
        "\tLog-Loss:\t", signif(object$log_loss, digits), "\n",
        "\tType of Posterior: ", object$type, "\n",
        "\tConfusion Matrix:", "\n")
    print(object$contingency_table)

    cat("\nPosterior Estimated Params: \n")

    print(object$bpwpm_params$betas)
    cat("\n")
    print(object$bpwpm_params$w)

    if(verb){
        cat("\nNodes\n")
        print(object$bpwpm_params$tau)

        # cat("\nF_Transformation\n")
        # print(object$bpwpm_params$estimated_F)
    }
}



