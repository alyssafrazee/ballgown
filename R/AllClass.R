#' Ballgown
#'
#' S4 class for storing and manipulating expression data from assembled transcriptomes
#'
#' @aliases Ballgown
#' 
#' @slot data tables containing expression data for genomic features (introns, exons, transcripts)
#' @slot structure genomic locations of features and their relationships to one another
#' @slot indexes tables connecting components of the assembly and providing other experimental 
#' information (e.g., phenotype data and locations of read alignment files)
#' @slot dirs directories holding data created by \code{tablemaker}
#' @slot mergedDate date the ballgown object was created
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
        data = "list",             # coverage data
        indexes = "list",          # reference information
        structure = "list",        # assembly information
        dirs = "character",        # directories where ballgown data is stored
        mergedDate = "character"   # date the object was created
    )
)
