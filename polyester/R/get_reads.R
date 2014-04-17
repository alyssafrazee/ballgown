#' get sequencing reads from fragments
#'
#' simulate the sequencing process by returning the sequence of one or both ends of provided 
#' fragments
#' @param tFrags DNAStringSet representing fragments
#' @param readlen Read length.
#' @param paired If \code{FALSE}, return only the first \code{readlen} bases of each element of 
#' \code{tFrags} in the result; if \code{TRUE}, also return last \code{readlen} bases.
#' @export
#' @return DNAStringSet representing simulated RNA-seq reads
#' @seealso \code{\link{simulate_experiment}}, \code{\link{simulate_experiment_countmat}}
get_reads = function(tFrags, readlen, paired = TRUE){
  
    # when fragments are shorter than reads:
    isShort = (width(tFrags) <= readlen)
    isLong = !isShort
      
    if(paired) {
      
        if(sum(isShort) > 0){
            x = tFrags[isShort]
            names(x) = paste0(seq(along=x), "a")
            rc = reverseComplement(x)
            names(rc) = paste0(seq(along=x), "b")
            out = c(x,rc)
            outShort = out[order(names(out))] # puts pairs of reads next to each other
            names(outShort) = paste0(rep(names(tFrags)[isShort],each=2))
        }
    
        if(sum(isLong) > 0){
            x = tFrags[isLong]
            lr = subseq(x, start=1, end=readlen)
            names(lr) = paste0(seq(along=x), "a")
            rr = subseq(x, start=(width(x)-readlen+1), end=width(x))
            rr = reverseComplement(rr)
            names(rr) = paste0(seq(along=x), "b")
            out = c(lr, rr)
            outLong = out[order(names(out))] # puts pairs of reads next to each other
            names(outLong) = paste0(rep(names(tFrags)[isLong],each=2))    
        }
      
        if(sum(isLong) > 0 & sum(isShort)) {
            theReads = c(outLong, outShort)
        } else if(sum(isLong) > 0) {
            theReads = outLong
        } else {
            theReads = outShort
        }
      
        return(theReads)
    
    } else { #  single end
      
      theReads = tFrags
      theReads[isLong] = subseq(tFrags[isLong], start=1, end=readlen)
      return(theReads)
    }
    
}
