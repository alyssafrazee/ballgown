#' read in GTF/GFF file as a data frame
#'
#' @param gffFile name of GTF/GFF on disk
#' @param nrows optional argument; number of rows to read in (default -1, which means read all rows)
#' @return data frame representing the GTF/GFF file
#' @seealso \code{\link{getAttributeField}} to extract data from "attributes" column; \link{http://useast.ensembl.org/info/website/upload/gff.html} for more information on GTF/GFF files.
#' @author Kasper Hansen
#' @export
gffRead <- function (gffFile, nrows = -1, verbose=FALSE) 
{
    if(verbose){
        cat("Reading ", gffFile, ": ", sep = "")    
    }
    gff = read.table(gffFile, sep = "\t", as.is = TRUE, quote = "", 
        header = FALSE, comment.char = "#", nrows = nrows, colClasses = c("character", 
            "character", "character", "integer", "integer", "character", 
            "character", "character", "character"))
    colnames(gff) = c("seqname", "source", "feature", "start", 
        "end", "score", "strand", "frame", "attributes")
    if(verbose){
        cat("found", nrow(gff), "rows with classes:", paste(sapply(gff, 
            class), collapse = ", "), "\n")
    }
    stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
    return(gff)
}

# source: https://stat.ethz.ch/pipermail/bioconductor/2008-October/024669.html