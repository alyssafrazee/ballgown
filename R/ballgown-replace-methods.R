setReplaceMethod("indexes", "ballgown", function(x, value) {x@indexes <- value; x})
setReplaceMethod("data", "ballgown", function(x, value) {x@data <- value; x})

# not called ReplaceMethod, but works for defining pData <- 
setMethod("pData", "ballgown", function(x){
  return(indexes(x)$pData)
})
setMethod("pData<-", "ballgown", function(x, value) {x@indexes$pData <- value; x})
