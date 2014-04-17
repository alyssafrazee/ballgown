#' reverse-complement some fragments
#'
#' randomly reverse-complement half of the sequences in a DNAStringSet
#' @param tObj DNAStringSet representing sequences.
#' @param seed optional seed to set before randomly selecting the sequences to be 
#' reverse-complemented.
#' @export
#' @return DNAStringSet that is the same as \code{tObj}, but with about half the sequences 
#' reverse-complemented.
reverse_complement = function(tObj, seed=NULL){
    if(!is.null(seed)) set.seed(seed)
    strand = sample(c(0,1), length(tObj), replace=TRUE)
    tObj[strand==0] = reverseComplement(tObj[strand==0])
    return(tObj)
}


