# function to draw reads from negative binomial

NB = function(inds, basemeans, dispersion_param, seed=NULL){
    if(!is.null(seed)) set.seed(seed)
    numreads = sapply(inds, function(x){
        rnbinom(n = 1, mu = basemeans[x], size = dispersion_param)
    }, USE.NAMES=FALSE)
    numreads[numreads == 0] = 1
    return(numreads)
}
