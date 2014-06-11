#' get names of samples in a ballgown objects
#' 
#' @name sampleNames
#' @exportMethod sampleNames
#' @docType methods
#' @rdname sampleNames
#' @aliases sampleNames,ballgown-method
#' @param x a ballgown object
#' @return vector of sample IDs for \code{x}. If \code{pData} exists, samples in its rows correspond 
#'  to samples in \code{sampleNames(x)} (in order).
setMethod("sampleNames", "ballgown", function(x){
    return(names(x@dirs))
})
