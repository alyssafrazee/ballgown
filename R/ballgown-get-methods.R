# slot getters

#' extract structure components from ballgown objects
#' 
#' @name structure
#' @exportMethod structure
#' @docType methods
#' @rdname structure
#' @aliases structure,ballgown-method
#' @param x a ballgown object
#' @return list containing elements \code{intron}, \code{exon}, and \code{trans}. \code{exon} and 
#'  \code{intron} are \code{GRanges} objects, where each range is an exon or intron, and 
#'  \code{trans} is a \code{GRangesList} object, where each \code{GRanges} element is a set of exons
#'  representing a transcript.
setMethod("structure", "ballgown", function(x) x@structure)

#' extract expression components from ballgown objects
#' 
#' @name expr
#' @exportMethod expr
#' @docType methods
#' @rdname expr
#' @aliases expr,ballgown-method
#' @param x a ballgown object
#' @return list containing elements \code{intron}, \code{exon}, and \code{trans}, 
#'  which are feature-by-sample data frames of expression data.
setMethod("expr", "ballgown", function(x) x@expr)

#' extract the indexes from ballgown objects
#' 
#' @name indexes
#' @exportMethod indexes
#' @docType methods
#' @rdname indexes
#' @aliases indexes,ballgown-method
#' @param x a ballgown object
#' @return list containing elements \code{e2t}, \code{i2t}, \code{t2g}, 
#'  \code{bamfiles}, and \code{pData}, where \code{e2t} and \code{i2t} are data frames linking exons 
#'  and introns (respectively) to transcripts, \code{t2g} is a data frame linking transcripts to 
#'  genes, and \code{bamfiles} and \code{pData} are described in \code{?ballgown}.
setMethod("indexes", "ballgown", function(x) x@indexes)

#' extract paths to tablemaker output 
#' 
#' @name dirs
#' @exportMethod dirs
#' @docType methods
#' @rdname dirs
#' @aliases dirs,ballgown-method
#' @param x a ballgown object
setMethod("dirs", "ballgown", function(x) x@dirs)

#' extract package version & creation date from ballgown object
#' 
#' @name mergedDate
#' @exportMethod mergedDate
#' @docType methods
#' @rdname mergedDate
#' @aliases mergedDate,ballgown-method
#' @param x a ballgown object
setMethod("mergedDate", "ballgown", function(x) x@mergedDate)

#' extract phenotype data from a ballgown object
#' 
#' @name pData
#' @exportMethod pData
#' @docType methods
#' @rdname pData
#' @aliases pData,ballgown-method
#' @param x a ballgown object
#' @return sample-by-phenotype data frame
setMethod("pData", "ballgown", function(x){
    return(indexes(x)$pData)
})

