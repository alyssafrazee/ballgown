### define generics for ballgown class 

#' @name structure
#' @export
#' @docType methods
#' @rdname structure
setGeneric("structure", function(x) standardGeneric("structure"))

#' @name expr
#' @export
#' @docType methods
#' @rdname expr
setGeneric("expr", function(x) standardGeneric("expr"))

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
#' @importFrom Biobase sampleNames
setGeneric("sampleNames")

#' @name mergedDate
#' @export
#' @docType methods
#' @rdname mergedDate
setGeneric("mergedDate", function(x) standardGeneric("mergedDate"))

#' @name indexes<-
#' @export
#' @docType methods
#' @rdname indexes-replace
setGeneric("indexes<-", function(x, value) standardGeneric("indexes<-"))

#' @name expr<-
#' @export
#' @docType methods
#' @rdname expr-replace
setGeneric("expr<-", function(x, value) standardGeneric("expr<-"))

#' @name subset
#' @export
#' @docType methods
#' @rdname subset
#' @param ... further arguments to generic subset
setGeneric("subset", function(x, ...) standardGeneric("subset"))

#' @name pData
#' @export
#' @docType methods
#' @rdname pData
#' @importFrom Biobase pData
setGeneric("pData")

#' @name pData<-
#' @export
#' @docType methods
#' @rdname pData-replace
#' @importFrom Biobase pData<-
setGeneric("pData<-")

#' @name texpr
#' @export
#' @docType methods
#' @rdname texpr
setGeneric("texpr", function(x, meas='FPKM') standardGeneric("texpr"))

#' @name eexpr
#' @export
#' @docType methods
#' @rdname eexpr
setGeneric("eexpr", function(x, meas='rcount') standardGeneric("eexpr"))

#' @name iexpr
#' @export
#' @docType methods
#' @rdname iexpr
setGeneric("iexpr", function(x, meas='rcount') standardGeneric("iexpr"))

#' @name gexpr
#' @export
#' @docType methods
#' @rdname gexpr
setGeneric("gexpr", function(x) standardGeneric("gexpr"))

#' @name seqnames
#' @export
#' @docType methods
#' @rdname seqnames
#' @importFrom GenomeInfoDb seqnames
setGeneric("seqnames")

#' @name transcriptIDs
#' @export
#' @docType methods
#' @rdname transcriptIDs
setGeneric('transcriptIDs', function(x) standardGeneric('transcriptIDs'))

#' @name transcriptNames
#' @export
#' @docType methods
#' @rdname transcriptNames
setGeneric('transcriptNames', function(x) standardGeneric('transcriptNames'))

#' @name geneIDs
#' @export
#' @docType methods
#' @rdname geneIDs
setGeneric('geneIDs', function(x) standardGeneric('geneIDs'))

#' @name geneNames
#' @export
#' @docType methods
#' @rdname geneNames
setGeneric('geneNames', function(x) standardGeneric('geneNames'))



