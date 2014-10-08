# slot getters

#' extract structure components from ballgown objects
#' 
#' @name structure
#' @exportMethod structure
#' @docType methods
#' @rdname structure
#' @aliases structure,ballgown-method
#' @param x a ballgown object
#' @return list containing elements \code{intron}, \code{exon}, and 
#'  \code{trans}. \code{exon} and \code{intron} are \code{GRanges} objects,
#'  where each range is an exon or intron, and \code{trans} is a
#'  \code{GRangesList} object, where each \code{GRanges} element is a set of 
#'  exons representing a transcript.
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
#' @return list containing elements \code{intron}, \code{exon}, and 
#'  \code{trans}, which are feature-by-sample data frames of expression data.
#' 
#' @examples
#' data(bg)
#' names(expr(bg))
#' class(expr(bg))
#' dim(expr(bg)$exon)
#' 
#' @seealso \code{\link{texpr}}, \code{\link{gexpr}}, \code{\link{eexpr}}, 
#'   \code{\link{iexpr}}
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
#'  \code{bamfiles}, and \code{pData}, where \code{e2t} and \code{i2t} are data
#'  frames linking exons and introns (respectively) to transcripts, \code{t2g}
#'  is a data frame linking transcripts to genes, and \code{bamfiles} and 
#'  \code{pData} are described in \code{?ballgown}.
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
#' @return vector of sample IDs for \code{x}. If \code{pData} exists, samples in
#' its rows correspond to samples in \code{sampleNames(x)} (in order).
#' 
#' @examples
#' data(bg)
#' sampleNames(bg)
setMethod("sampleNames", "ballgown", function(object){
    return(names(object@dirs))
})

#' get numeric transcript IDs from a ballgown object
#'
#' @name transcriptIDs
#' @exportMethod transcriptIDs
#' @docType methods
#' @rdname transcriptIDs
#' @aliases transcriptIDs,ballgown-method
#' @param x a ballgown object
#' @return vector of numeric transcript IDs included in the ballgown object
#' 
#' @examples
#' data(bg)
#' transcriptIDs(bg)
setMethod("transcriptIDs", "ballgown", function(x){
    return(texpr(x, 'all')$t_id)
})

#' get transcript names from a ballgown object
#'
#' @name transcriptNames
#' @exportMethod transcriptNames
#' @docType methods
#' @rdname transcriptNames
#' @aliases transcriptNames,ballgown-method
#' @param x a ballgown object
#' @return vector of transcript names included in the ballgown object. If object
#' was created using Cufflinks/Tablemaker, these transcript names will be of the
#' form "TCONS_*". Return vector is named and ordered by corresponding numeric
#' transcript ID.
#' 
#' @examples
#' data(bg)
#' transcriptNames(bg)
setMethod("transcriptNames", "ballgown", function(x){
    ret = texpr(x, 'all')$t_name
    names(ret) = texpr(x, 'all')$t_id
    ret
})

#' get gene IDs from a ballgown object
#'
#' @name geneIDs
#' @exportMethod geneIDs
#' @docType methods
#' @rdname geneIDs
#' @aliases geneIDs,ballgown-method
#' @param x a ballgown object
#' @return named vector of gene IDs included in the ballgown object. If object 
#' was created using Tablemaker, these gene IDs will be of the form "XLOC_*".
#' Vector is named and ordered by corresponding numeric transcript ID.
#'
#' @details
#' This vector differs from that produced by geneNames in that geneIDs produces
#' names of loci created during the assembly process, not necessarily
#' annotated genes. 
#' @seealso \code{\link{geneNames}}
#' @examples
#' data(bg)
#' geneIDs(bg)
setMethod("geneIDs", "ballgown", function(x){
    ret = texpr(x, 'all')$gene_id
    names(ret) = texpr(x, 'all')$t_id
    ret
})

#' get gene names from a ballgown object
#'
#' @name geneNames
#' @exportMethod geneNames
#' @docType methods
#' @rdname geneNames
#' @aliases geneNames,ballgown-method
#' @param x a ballgown object
#' @return named vector of gene names included in the ballgown object, named
#' and ordered by corresponding numeric transcript ID. 
#'
#' @details
#' This vector differs from that produced by geneIDs in that geneNames produces
#' *annotated* gene names that correspond to assembled transcripts. The
#' return will be empty/blank/NA if the transcriptome assembly is de novo 
#' (i.e., was not compared to an annotation before the ballgown object was
#' created). See \code{\link{getGenes}} for matching transcripts to gene names.
#' Some entries of this vector will be empty/blank/NA if the corresponding 
#' transcript did not overlap any annotated genes. 
#' 
#' @examples
#' data(bg)
#' # this is a de novo assembly, so it does not contain gene info as it stands
#' # but we can add it:
#' annot = system.file('extdata', 'annot.gtf.gz', package='ballgown')
#' gnames = getGenes(annot, structure(bg)$trans, UCSC=FALSE)
#' gnames_first = lapply(gnames, function(x) x[1]) #just take 1 overlapping gene
#' expr(bg)$trans$gene_name = gnames_first
#' 
#' # now we can extract these gene names:
#' geneNames(bg)
#' 
#' @seealso \code{\link{geneIDs}}
setMethod("geneNames", "ballgown", function(x){
    ret = texpr(x, 'all')$gene_name
    names(ret) = texpr(x, 'all')$t_id
    ret
})

#' get sequence (chromosome) names from ballgown object
#'
#' @name seqnames
#' @exportMethod seqnames
#' @docType methods
#' @rdname seqnames
#' @aliases seqnames,ballgown-method
#' @param x a ballgown object
#' @return vector of sequence (i.e., chromosome) names included in the 
#' ballgown object
#' 
#' @examples
#' data(bg)
#' seqnames(bg)
setMethod("seqnames", "ballgown", function(x){
    return(unique(texpr(x,'all')$chr))
})