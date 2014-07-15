#' nice wrapper for string splitting and subsetting
#'
#' @param x character vector whose elements are to be split
#' @param pattern what to split by
#' @param slot which slot from each element of x you'd like returned (default 1)
#' @param ... other arguments to strsplit
#' @details internal use only 
#' @author Andrew Jaffe
#' @examples
#' x = c('sample1.FPKM', 'sample1.cov', 'sample2.FPKM', 'sample2.cov')
#' ss(x, pattern='\\.', slot=1) # sample IDs from x
ss = function(x, pattern, slot=1, ...){
    sapply(strsplit(x,pattern,...), "[", slot)  
} 
