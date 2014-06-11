### define generics for ballgown class 

#' @name structure
#' @export
#' @docType methods
#' @rdname structure
setGeneric("structure", function(x) standardGeneric("structure"))

#' @name data
#' @export
#' @docType methods
#' @rdname data
setGeneric("data", function(x) standardGeneric("data"))

#' @name indexes
#' @export
#' @docType methods
#' @rdname indexes
setGeneric("indexes", function(x) standardGeneric("indexes"))

#' @name dirs
#' @export
#' @docType methods
#' @rdname dirs
setGeneric("dirs", function(x) standardGeneric("dirs"))

#' @name sampleNames
#' @export
#' @docType methods
#' @rdname sampleNames
setGeneric("sampleNames", function(x) standardGeneric("sampleNames"))

#' @name mergedDate
#' @export
#' @docType methods
#' @rdname mergedDate
setGeneric("mergedDate", function(x) standardGeneric("mergedDate"))

#' @name indexes<-
#' @export
#' @docType methods
#' @rdname ballgown-replace-methods
setGeneric("indexes<-", function(x, value) standardGeneric("indexes<-"))

#' @name data<-
#' @export
#' @docType methods
#' @rdname ballgown-replace-methods
setGeneric("data<-", function(x, value) standardGeneric("data<-"))

#' @name subset
#' @export
#' @docType methods
#' @rdname subset
setGeneric("subset", function(x, ...) standardGeneric("subset"))

#' @name pData
#' @export
#' @docType methods
#' @rdname pData
setGeneric("pData", function(x) standardGeneric("pData"))

#' @name pData<-
#' @exportMethod pData<-
#' @docType methods
#' @rdname ballgown-replace-methods
setGeneric("pData<-", function(x, value) standardGeneric("pData<-"))

#' @name texpr
#' @export
#' @docType methods
#' @rdname texpr
setGeneric("texpr", function(x, ...) standardGeneric("texpr"))

#' @name eexpr
#' @export
#' @docType methods
#' @rdname eexpr
setGeneric("eexpr", function(x, ...) standardGeneric("eexpr"))

#' @name iexpr
#' @exportMethod iexpr
#' @docType methods
#' @rdname iexpr
setGeneric("iexpr", function(x, ...) standardGeneric("iexpr"))

#' @name gexpr
#' @export
#' @docType methods
#' @rdname gexpr
setGeneric("gexpr", function(x) standardGeneric("gexpr"))




