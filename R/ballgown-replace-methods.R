#' Replace method for indexes slot in ballgown objects
#' 
#' @name indexes<-
#' @exportMethod indexes<-
#' @docType methods
#' @aliases indexes<-,ballgown-method
#' @rdname indexes-replace
#' @param x a ballgown object
#' @param value the updated value for \code{indexes(x)} or a subcomponent
setReplaceMethod("indexes", "ballgown", function(x, value) {x@indexes <- value; x})

#' Replacement method for expr slot in ballgown objects
#' 
#' @name expr<-
#' @exportMethod expr<-
#' @docType methods
#' @aliases expr<-,ballgown-method
#' @rdname expr-replace
#' @param x a ballgown object
#' @param value the updated value for \code{expr(x)} or a subcomponent
setReplaceMethod("expr", "ballgown", function(x, value) {x@expr <- value; x})

# setReplaceMethod did not work (since it's a subcomponent of indexes?)

#' Replacement method for pData slot in ballgown objects
#' 
#' @name pData<-
#' @exportMethod pData<-
#' @docType methods
#' @aliases pData<-,ballgown,ANY-method
#' @rdname pData-replace
#' @param object a ballgown object
#' @param value the updated value for \code{pData(x)}.
setMethod("pData<-", "ballgown", function(object, value) {
    object@indexes$pData <- value
    ### TODO: check validity of provided pData
    object}
)
