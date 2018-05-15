log_loss <- function(y, p, eps = 1e-15, verb = TRUE){
    "Function that evaluates if the log loss of a binary probability model is withina a certain threshold
    
    INPUTS: 
    y:= Vector of observed binary outcomes (integer - n)
    p:= Vector of fitted probabilities ie: p = P(y = 1) (numeric - n)
    eps:= Error for evaluating the logarithm function. Important for numerical stability (numeric)
    verb:= Boolean for printing information (boolen)
    
    OUTPUT:
    LogLoss:= value of the LogLoss function for binary outcomes (numeric)
    "
    p_corrected = pmin(pmax(p, eps), 1-eps)
    ll <- - sum (y * log(p_corrected) + (1 - y) * log(1 - p_corrected))/length(y)
    if(verb){
        cat("\nLog_Loss: ",ll, sep ="")
    }
    return(- sum (y * log(p_corrected) + (1 - y) * log(1 - p_corrected))/length(y))
}


