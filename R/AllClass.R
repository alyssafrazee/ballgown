#' Ballgown
#'
#' S4 class for storing and manipulating expression data from assembled 
#' transcriptomes
#'
#' @aliases Ballgown
#' 
#' @slot expr tables containing expression data for genomic features (introns, 
#' exons, transcripts)
#' @slot structure genomic locations of features and their relationships to one 
#'   another
#' @slot indexes tables connecting components of the assembly and providing 
#'   other experimental 
#' information (e.g., phenotype data and locations of read alignment files)
#' @slot dirs directories holding data created by \code{tablemaker}
#' @slot mergedDate date the ballgown object was created
#' @slot meas which expression measurement(s) the object contains in its data 
#'   slot. Vector of one or more of "rcount", "ucount", "mrcount", "cov", 
#'   "cov_sd", "mcov", "mcov_sd", or "FPKM", if Tablemaker output is used, or 
#'   one of "TPM" or "FPKM" if RSEM output is used. Can also be "all" for all 
#'   measurements. See vignette for details.
#' @slot RSEM TRUE if object was made from RSEM output, FALSE if object was made
#'   from Tablemaker/Cufflinks output.
#' 
#' @name ballgown-class
#' 
#' @rdname ballgown-class
#' 
#' @exportClass ballgown
#' 
#' @author Alyssa Frazee, Leonardo Collado-Torres, Jeff Leek
#' @examples
#'   data(bg)
#'   class(bg) #"ballgown"
#'   dim(bg@@expr$exon)
#'   bg@@structure$exon
#'   head(bg@@indexes$t2g)
#'   head(bg@@dirs)
#'   bg@@mergedDate
#'   bg@@meas
#'   bg@@RSEM
setClass("ballgown", 
    representation(
        expr = "list",             # coverage data
        indexes = "list",          # reference information
        structure = "list",        # assembly information
        dirs = "character",        # directories where ballgown data is stored
        mergedDate = "character",  # date the object was created
        meas = "character",        # what measurements are in the object
        RSEM = "logical"           # TRUE if made from RSEM output
    )
)
