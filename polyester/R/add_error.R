### function to add error to reads/fragments

# Q: should this go after we actually select reads? will that matter? (for speed and/or sensibleness?)
#'<brief desc>
#'
#'<full description>
#' @param tFrags <what param does>
#' @param  error_rate = 0.005 <what param does>
#' @export
#' @keywords
#' @seealso
#' @return
#' @alias
#' @examples \dontrun{
#'
#'}
add_error = function(tFrags, error_rate = 0.005){
    adj_error = error_rate*4/3 # based on random read selection
    
    allSeq = unlist(tFrags)
    insertLocs = Rle(sample(c(TRUE,FALSE), size = length(allSeq), 
           replace=TRUE, prob = c(adj_error, 1-adj_error)))
  
    newletters = DNAString(paste(sample(c("A", "C", "G", "T"), sum(insertLocs), replace=TRUE), collapse="") )
    allSeq = replaceLetterAt(allSeq, insertLocs, newletters)
    
    eFrags = DNAStringSet(allSeq, start=c(1, cumsum(width(tFrags))[-length(tFrags)]), width=width(tFrags))
    names(eFrags) = names(tFrags)
    return(eFrags)
}
 
