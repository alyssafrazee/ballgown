#'  Draw nonzero negative binomial random numbers
#'
#' @param  basemeans vector of means, one per draw
#' @param  dispersion_param vector of dispersion (\code{size}) parameters, one per draw
#' @param  seed optional seed to set before drawing
NB = function(basemeans, dispersion_param, seed=NULL){
    if(!is.null(seed)) set.seed(seed)
    numreads = rnbinom(n = length(basemeans), mu = basemeans, size = dispersion_param)
    numreads[numreads == 0] = 1
    return(numreads)
}
