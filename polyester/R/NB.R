#'  Negative binomial draws
#'
#'<full description>
#' @param  basemeans <what param does>
#' @param  dispersion_param <what param does>
#' @param  seed=NULL <what param does>
#' @export
#' @examples \dontrun{
#'
#'}
NB = function(basemeans, dispersion_param, seed=NULL){
    if(!is.null(seed)) set.seed(seed)
    numreads = rnbinom(n = length(basemeans), mu = basemeans, size = dispersion_param)
    numreads[numreads == 0] = 1
    return(numreads)
}
