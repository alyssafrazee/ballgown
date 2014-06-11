#' extract a specific field of the "attributes" column of a data frame created from a GTF/GFF file
#'
#' @param x vector representing the "attributes" column of GTF/GFF file
#' @param field name of the field you want to extract from the "attributes" column
#' @param attrsep separator for the fields in the attributes column.  Defaults to '; ', the 
#' separator for GTF files outputted by Cufflinks.
#' 
#' @return vector of nucleotide positions included in the transcript
#' 
#' @seealso \code{\link{gffRead}} for creating a data frame from a GTF/GFF file, and 
#' \url{http://useast.ensembl.org/info/website/upload/gff.html} for specifics of the GFF/GTF file 
#' format.
#' 
#' @author Wolfgang Huber, in \code{davidTiling}
#' 
#' @export
#' @examples \dontrun{
#' ## not run:
#' gffdata = gffRead('transcripts.gtf')
#' transcriptID = getAttributeField(gffdata$V9, field = "transcript_id")
#' gffdata$transcriptID = transcriptID
#' }
getAttributeField <- function (x, field, attrsep = "; ") 
{
    s = strsplit(x, split = attrsep, fixed = TRUE)
    sapply(s, function(atts) {
        a = strsplit(atts, split = " ", fixed = TRUE)
        m = match(field, sapply(a, "[", 1))
        if (!is.na(m)) {
            rv = a[[m]][2]
        }
        else {
            rv = as.character(NA)
        }
        return(rv)
    })
}

### needs attribution
### https://stat.ethz.ch/pipermail/bioconductor/2008-October/024669.html
### this function comes from the davidTiling package, but I changed the default attribute separator
### don't want to import the whole package, but would like to include this function in ballgown