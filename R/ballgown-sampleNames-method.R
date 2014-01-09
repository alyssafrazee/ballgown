setMethod("sampleNames", "ballgown", function(x){
    return(names(x@dirs))
})
