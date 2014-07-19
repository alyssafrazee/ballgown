#' Replace method for indexes slot in ballgown objects
#' 
#' @name indexes<-
#' @exportMethod indexes<-
#' @docType methods
#' @aliases indexes<-,ballgown-method
#' @rdname indexes-replace
#' @param x a ballgown object
#' @param value the updated value for \code{indexes(x)} or a subcomponent
#' 
#' @examples
#' data(bg)
#' indexes(bg)$bamfiles = paste0('/path/to/bamfolder/', 
#'   sampleNames(bg), '_accepted_hits.bam')
setReplaceMethod("indexes", "ballgown", function(x, value){
    x@indexes <- value; x})

#' Replacement method for expr slot in ballgown objects
#' 
#' @name expr<-
#' @exportMethod expr<-
#' @docType methods
#' @aliases expr<-,ballgown-method
#' @rdname expr-replace
#' @param x a ballgown object
#' @param value the updated value for \code{expr(x)} or a subcomponent
#' 
#' @examples
#' data(bg)
#' n = ncol(bg@@expr$trans)
#' #multiply all transcript expression measurements by 10:
#' bg@@expr$trans[,11:n] = 10*bg@@expr$trans[11:n] 
setReplaceMethod("expr", "ballgown", function(x, value) {x@expr <- value; x})


#' Replacement method for pData slot in ballgown objects
#' 
#' @name pData<-
#' @exportMethod pData<-
#' @docType methods
#' @aliases pData<-,ballgown,ANY-method
#' @rdname pData-replace
#' @param object a ballgown object
#' @param value the updated value for \code{pData(x)}.
#' 
#' @examples
#' # add "timepoint" covariate to ballgown object:
#' data(bg) # already contains pData
#' pData(bg) = data.frame(pData(bg), timepoint=rep(1:10, 2))
#' head(pData(bg))
setMethod("pData<-", "ballgown", function(object, value) {
    object@indexes$pData <- value
    object}
)
