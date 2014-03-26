# generate fragments from DNAStringSetList of transcripts
# new feature: variable fragment length

generate_fragments = function(tObj, fraglen, fragsd = 25){
    L = width(tObj)
    fraglens = round(rnorm(L, mean=fraglen, sd=fragsd)) #add variable fragment lengths
    s = which(fraglens < L)
    tObj[s] = subseq(tObj[s],start = floor(runif(length(s), min = rep(1,length(s)), 
                                                max = L[s]-fraglens[s])), width=fraglens[s])
    return(tObj)
}

### probably don't want all fragments of the same length for same transcript


