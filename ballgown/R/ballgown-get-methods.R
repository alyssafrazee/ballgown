# slot getters

setMethod("structure", "ballgown", function(x) x@structure)
setMethod("data", "ballgown", function(x) x@data)
setMethod("indexes", "ballgown", function(x) x@indexes)
setMethod("dirs", "ballgown", function(x) x@dirs)
setMethod("mergedDate", "ballgown", function(x) x@mergedDate)
