#' @importFrom greybox graphmaker
#' @export
plot.counter <- function(x, ...){
    if(any(!is.na(c(x$lower,x$upper)))){
        graphmaker(x$actuals, x$forecast, x$fitted, x$lower, x$upper, x$level, ...);
    }
    else{
        graphmaker(x$actuals, x$forecast, x$fitted, ...);
    }
}

#' @export
print.counter <- function(x, ...){
    cat("The model constructed: "); cat(x$model);
    if(any(!is.na(x$holdout))){
        cat("\nAccuracy:\n");
        print(x$accuracy[c("sMAE","sMSE","sCE","RelMAE","RelMSE","RelAME","sPIS")]);
    }
}

#' @importFrom smooth pls
#' @importFrom stats dnbinom dpois window
#' @export
pls.counter <- function(object, holdout=NULL, ...){
    if(is.null(holdout)){
        if(is.na(object$holdout)){
            stop("We need the values from the holdout in order to proceed...", call.=FALSE);
        }
        else{
            holdout <- object$holdout;
        }
    }
    
    if(object$model=="HSP"){
        return(sum(dpois(x=holdout, lambda=object$forecast, log=TRUE)));
    }
    else if(object$model=="NegBin"){
        negBinSize <- window(object$dispersion, start=start(object$forecast));
        return(sum(dnbinom(x=holdout, size=negBinSize, mu=object$forecast, log=TRUE)));
    }
}