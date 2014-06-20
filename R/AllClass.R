#' Ballgown
#'
#' S4 class for storing and manipulating expression data from assembled transcriptomes
#'
#' @aliases Ballgown
#' 
#' @slot expr tables containing expression data for genomic features (introns, exons, transcripts)
#' @slot structure genomic locations of features and their relationships to one another
#' @slot indexes tables connecting components of the assembly and providing other experimental 
#' information (e.g., phenotype data and locations of read alignment files)
#' @slot dirs directories holding data created by \code{tablemaker}
#' @slot mergedDate date the ballgown object was created
#' @slot meas which expression measurement(s) the object contains in its data slot. Vector of one or
#'   more of "rcount", "ucount", "mrcount", "cov", "cov_sd", "mcov", "mcov_sd", or "FPKM". See 
#'   vignette for details.
#' 
#' @name ballgown-class
#' 
#' @rdname ballgown-class
#' 
#' @exportClass ballgown
#' 
#' @author Alyssa Frazee, Leonardo Collado Torres, Jeff Leek
setClass("ballgown", 
    representation(
        expr = "list",             # coverage data
        indexes = "list",          # reference information
        structure = "list",        # assembly information
        dirs = "character",        # directories where ballgown data is stored
        mergedDate = "character",  # date the object was created
        meas = "character"         # flags object as only containing certain measurements (eg FPKM)
    )
)
