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
#'
#' @examples
#' data(bg)
#' names(structure(bg))
#' class(structure(bg))
#' structure(bg)$exon
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
#' 
#' @examples
#' data(bg)
#' names(expr(bg))
#' class(expr(bg))
#' dim(expr(bg)$exon)
#' 
#' @seealso \code{\link{texpr}}, \code{\link{gexpr}}, \code{\link{eexpr}}, \code{\link{iexpr}}
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
#' 
#' @examples
#' data(bg)
#' names(indexes(bg))
#' class(indexes(bg))
#' head(indexes(bg)$t2g)
setMethod("indexes", "ballgown", function(x) x@indexes)

#' extract paths to tablemaker output 
#' 
#' @name dirs
#' @exportMethod dirs
#' @docType methods
#' @rdname dirs
#' @aliases dirs,ballgown-method
#' @param x a ballgown object
#' 
#' @examples
#' data(bg)
#' dirs(bg)
setMethod("dirs", "ballgown", function(x) x@dirs)

#' extract package version & creation date from ballgown object
#' 
#' @name mergedDate
#' @exportMethod mergedDate
#' @docType methods
#' @rdname mergedDate
#' @aliases mergedDate,ballgown-method
#' @param x a ballgown object
#' 
#' @examples
#' data(bg)
#' mergedDate(bg)
setMethod("mergedDate", "ballgown", function(x) x@mergedDate)

#' extract phenotype data from a ballgown object
#' 
#' @name pData
#' @exportMethod pData
#' @docType methods
#' @rdname pData
#' @aliases pData,ballgown-method
#' @param object a ballgown object
#' @return sample-by-phenotype data frame
#' 
#' @examples
#' data(bg)
#' pData(bg)
setMethod("pData", "ballgown", function(object){
    return(indexes(object)$pData)
})

#' get names of samples in a ballgown objects
#' 
#' @name sampleNames
#' @exportMethod sampleNames
#' @docType methods
#' @rdname sampleNames
#' @aliases sampleNames,ballgown-method
#' @param object a ballgown object
#' @return vector of sample IDs for \code{x}. If \code{pData} exists, samples in its rows correspond 
#'  to samples in \code{sampleNames(x)} (in order).
#' 
#' @examples
#' data(bg)
#' sampleNames(bg)
setMethod("sampleNames", "ballgown", function(object){
    return(names(object@dirs))
})

