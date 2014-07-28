# the show method:
setMethod("show", "ballgown", 
    function(object)
        cat(class(object), "instance with", length(structure(object)$trans), 
            "assembled transcripts and", length(object@dirs), "samples\n")
)
