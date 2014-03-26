# function to get reads from list of fragments

get_reads = function(tFrags, readlen, paired = TRUE){
  
    # when fragments are shorter than reads:
    isShort = (width(tFrags) <= readlen)
    isLong = !isShort
      
    if(paired) {
      
      if(sum(isShort) > 0){
        x = tFrags[isShort]
        names(x) = paste0(seq(along=x), "a")
        rc = reverseComplement(x)
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
      theReads[isLong] = subseq(tFrags[isLong],start=1, end=readlen)
      return(theReads)
    }
    
}
