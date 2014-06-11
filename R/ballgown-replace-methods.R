#' @name indexes<-
#' @exportMethod indexes<-
#' @docType methods
#' @rdname ballgown-replace-methods
#' @param x a ballgown object
#' @param value the updated value for \code{indexes(x)} or a subcomponent
setReplaceMethod("indexes", "ballgown", function(x, value) {x@indexes <- value; x})

#' @name data<-
#' @exportMethod data<-
#' @docType methods
#' @rdname ballgown-replace-methods
#' @param x a ballgown object
#' @param value the updated value for \code{data(x)} or a subcomponent
setReplaceMethod("data", "ballgown", function(x, value) {x@data <- value; x})

# setReplaceMethod did not work (since it's a subcomponent of indexes?)
#' @name pData<-
#' @exportMethod pData<-
#' @docType methods
#' @aliases pData<-,ballgown-method
#' @rdname ballgown-replace-methods
#' @param x a ballgown object
#' @param value the updated value for \code{pData(x)}.
setMethod("pData<-", "ballgown", function(x, value) {
    x@indexes$pData <- value
    ### TODO: check validity of provided pData
    x}
)
