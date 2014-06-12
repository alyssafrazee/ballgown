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

#' Replacement method for data slot in ballgown objects
#' 
#' @name data<-
#' @exportMethod data<-
#' @docType methods
#' @aliases data<-,ballgown-method
#' @rdname data-replace
#' @param x a ballgown object
#' @param value the updated value for \code{data(x)} or a subcomponent
setReplaceMethod("data", "ballgown", function(x, value) {x@data <- value; x})

# setReplaceMethod did not work (since it's a subcomponent of indexes?)

#' Replacement method for pData slot in ballgown objects
#' 
#' @name pData<-
#' @exportMethod pData<-
#' @docType methods
#' @aliases pData<-,ballgown-method
#' @rdname pData-replace
#' @param x a ballgown object
#' @param value the updated value for \code{pData(x)}.
setMethod("pData<-", "ballgown", function(x, value) {
    x@indexes$pData <- value
    ### TODO: check validity of provided pData
    x}
)
