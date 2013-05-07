gffRead <- function (gffFile, nrows = -1) 
{
    cat("Reading ", gffFile, ": ", sep = "")
    gff = read.table(gffFile, sep = "\t", as.is = TRUE, quote = "", 
        header = FALSE, comment.char = "#", nrows = nrows, colClasses = c("character", 
            "character", "character", "integer", "integer", "character", 
            "character", "character", "character"))
    colnames(gff) = c("seqname", "source", "feature", "start", 
        "end", "score", "strand", "frame", "attributes")
    cat("found", nrow(gff), "rows with classes:", paste(sapply(gff, 
        class), collapse = ", "), "\n")
    stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
    return(gff)
}