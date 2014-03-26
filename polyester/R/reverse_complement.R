# reverse complement each fragment in a DNAStringsList (or whatever) of fragments

reverse_complement = function(tObj){
    strand = sample(c(0,1), length(tObj), replace=TRUE)
    tObj[strand==0] = reverseComplement(tObj[strand==0])
    return(tObj)
}


