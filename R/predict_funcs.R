predict.bpwpm <- function(object, new_Y, new_X,
                         thin = 0, burn_in = 0, type = 'mean', ...) {

    if(!('bpwpm' %in% class(object))){
        error("Object not of the class bpwpm")
        geterrmessage()
    }

    post_params <- posterior_params(object, thin, burn_in, type)
    p <- posterior_probs(new_X, post_params)

    model_predict <- list(info = object$info,
                          type = type,
                          params = post_params,
                          contingency_table = contingency_table(new_Y, p),
                          accuracy = accurracy(new_Y, p),
                          log_loss = log_loss(new_Y, p, verb = FALSE),
                          X = new_X,
                          Y = new_Y)

    class(model_predict) <- "bpwpm_prediction"

    return(model_predict)
}

