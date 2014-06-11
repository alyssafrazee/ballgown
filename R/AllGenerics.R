### define generics for ballgown class 

#' methods for objects of class \code{ballgown}
#'
#' The following methods are available for the S4 class \code{ballgown}.
#' @name ballgown-methods
#' @aliases structure
#' @exportMethod structure
#' @docType methods
#' @rdname ballgown-methods
#' @param x ballgown object
#' @return for \code{structure}: list containing elements \code{intron} (GRanges), \code{exon} 
#' (GRanges), and \code{trans} (GRangesList), denoting genomic positions of exons, introns, and 
#' transcripts (represented as sets of exons).
setGeneric("structure", function(x) standardGeneric("structure"))

#' @name data
#' @exportMethod data
#' @docType methods
#' @rdname ballgown-methods
#' @return for \code{data}: list containing elements \code{intron}, \code{exon}, and \code{trans} 
#' (all data frames) -- feature-by-sample expression tables.
setGeneric("data", function(x) standardGeneric("data"))

#' @name indexes
#' @exportMethod indexes
#' @docType methods
#' @rdname ballgown-methods
#' @return for \code{indexes}: list containing elements \code{e2t}, \code{i2t}, \code{t2g}, 
#' \code{bamfiles}, and \code{pData}, where \code{e2t} and \code{i2t} are data frames linking exons 
#' and introns (respectively) to transcripts, \code{t2g} is a data frame linking transcripts to 
#' genes, and \code{bamfiles} and \code{pData} are described at the \code{link{ballgown}} 
#' constructor help page.
setGeneric("indexes", function(x) standardGeneric("indexes"))

#' @name dirs
#' @exportMethod dirs
#' @docType methods
#' @rdname ballgown-methods
#' @return for \code{dirs}: paths to the on-disk directories holding the data (created with 
#' \code{tablemaker}) used to create \code{x}
setGeneric("dirs", function(x) standardGeneric("dirs"))

#' @name sampleNames
#' @exportMethod sampleNames
#' @docType methods
#' @rdname ballgown-methods
#' @return for \code{sampleNames}: names of samples included in the dataset.  Matches folder names 
#' of \code{dirs(x)}.
setGeneric("sampleNames", function(x) standardGeneric("sampleNames"))

#' @name mergedDate
#' @exportMethod mergedDate
#' @docType methods
#' @rdname ballgown-methods
#' @return for \code{mergedDate}: the date \code{x} was created
setGeneric("mergedDate", function(x) standardGeneric("mergedDate"))

#' @name indexes<-
#' @exportMethod indexes<-
#' @docType methods
#' @rdname ballgown-methods
#' @param value the updated value for a ballgown object component
setGeneric("indexes<-", function(x, value) standardGeneric("indexes<-"))

#' @name data<-
#' @exportMethod data<-
#' @docType methods
#' @rdname ballgown-methods
setGeneric("data<-", function(x, value) standardGeneric("data<-"))

#' @name subset
#' @exportMethod subset
#' @docType methods
#' @rdname ballgown-methods
#' @details To use \code{subset}, you must provide the \code{cond} argument as a string representing
#' a logical expression specifying your desired subset. The subset expression can either involve 
#' column names of \code{texpr(x, "all")} (if \code{genomesubset} is \code{TRUE}) or of 
#' \code{pData(x)} (if \code{genomesubset} is \code{FALSE}). For example, if you wanted a ballgown 
#' object for only chromosome 22, you might call \code{subset(x, "chr == 'chr22'")}. 
#' (Be sure to handle quotes within character strings appropriately). Note that \code{genomesubset} 
#' is \code{TRUE} by default.
setGeneric("subset", function(x, ...) standardGeneric("subset"))

#' @name pData
#' @exportMethod pData
#' @docType methods
#' @rdname ballgown-methods
setGeneric("pData", function(x) standardGeneric("pData"))

#' @name pData<-
#' @exportMethod pData<-
#' @docType methods
#' @rdname ballgown-methods
setGeneric("pData<-", function(x, value) standardGeneric("pData<-"))

#' @name texpr
#' @exportMethod texpr
#' @docType methods
#' @rdname ballgown-methods
#' @param ... for \code{subset}: arguments are \code{cond}, a string giving a subset condition 
#' (see details) and \code{genomesubsest}, which is \code{TRUE} if you want a ballgown object for 
#' only part of the genome, and \code{FALSE} if you want a ballgown object containing only some of 
#' the samples in the experiment. \code{genomesubset} is \code{TRUE} by default. For \code{*expr} 
#' methods: one of \code{'cov'}, \code{'FPKM'}, \code{'rcount'}, \code{'ucount'}, \code{'mrcount'},
#' \code{'cov_sd'}, \code{'mcov'}, or \code{'mcov_sd'}, depending on which type of expression 
#' measurement is desired.  Leave \code{...} blank to select \code{"FPKM"} for \code{texpr}, or
#' \code{"rcount"} for \code{eexpr} or \code{iexpr}.
setGeneric("texpr", function(x, ...) standardGeneric("texpr"))

#' @name eexpr
#' @export
#' @docType methods
#' @rdname eexpr
setGeneric("eexpr", function(x, ...) standardGeneric("eexpr"))

#' @name iexpr
#' @exportMethod iexpr
#' @docType methods
#' @rdname ballgown-methods
setGeneric("iexpr", function(x, ...) standardGeneric("iexpr"))

#' @name gexpr
#' @exportMethod gexpr
#' @docType methods
#' @rdname ballgown-methods
#' @param meas expression measurement to extract (for use with \code{*expr} methods).  Defaults to 
#' FPKM for \code{texpr} and \code{rcount} for \code{eexpr} and \code{iexpr}. Specifying 
#' \code{"all"} will return all expression measurements as well as extra feature-level information.
#' @return for \code{*expr} methods: a feature-by-sample table with the specified expression
#' measurement(s) in the cells. Return is a data frame if \code{meas} is \code{"all"} and a matrix 
#' otherwise.
setGeneric("gexpr", function(x) standardGeneric("gexpr"))




