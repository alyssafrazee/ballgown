### function to add error to reads/fragments

# Q: should this go after we actually select reads? will that matter? (for speed and/or sensibleness?)
add_error = function(tFrags, error_rate = 0.005){
    adj_error = error_rate*4/3 # based on random read selection
    
    allSeq = unlist(tFrags)
    insertLocs = Rle(sample(c(TRUE,FALSE), size = length(allSeq), 
           replace=TRUE, prob = c(adj_error, 1-adj_error)))
  
    newletters = DNAString(paste(sample(c("A", "C", "G", "T"), sum(insertLocs), replace=TRUE), collapse="") )
    allSeq = replaceLetterAt(allSeq, insertLocs, newletters)
    
    eFrags = DNAStringSet(allSeq, start=c(1, (cumsum(width(tFrags))+1)[-length(tFrags)]), width=width(tFrags))
    names(eFrags) = names(tFrags)
    return(eFrags)
}
 
