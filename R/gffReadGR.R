#' read in gtf file as GRanges object
#' 
#' (very) light wrapper for rtracklayer::import
#'  
#' @param gtf name of GTF/GFF file on disk
#' @param splitByTranscript if \code{TRUE}, return a \code{GRangesList} of transcripts; otherwise 
#' return a \code{GRanges} object containing all genomic features in \code{gtf}. Default 
#' \code{FALSE}.
#' @param identifier name of transcript identifier column of \code{attributes} field in \code{gtf}. 
#' Default \code{"transcript_id"}. Only used if \code{splitByTranscript} is \code{TRUE}.
#' @param sep field separator in the \code{attributes} field of \code{gtf}. Default \code{"; "} 
#' (semicolon + space). Only used if \code{splitByTranscript} is \code{TRUE}.
#' 
#' @return if \code{splitByTranscript} is \code{FALSE}, an object of class \code{GRanges} 
#' representing the genomic features in \code{gtf}. If \code{splitByTranscript} is TRUE, an object 
#' of class \code{GRangesList}, where each element is a \code{GRanges} object corresponding to an 
#' annotated transcript (designated in \code{names}).
#' 
#' @seealso \code{\link{gffRead}} for reading in a GTF file as a data frame rather than a 
#' \code{GRanges}/\code{GRangesList} object.
#' 
#' @author Alyssa Frazee
#' 
#' @export
gffReadGR = function(gtf, splitByTranscript=FALSE, identifier='transcript_id', sep='; '){
    con = file(gtf)
    ret = import(con, format='GFF')
    if(splitByTranscript){
        ret = ret[ret$type == 'exon', ]
        transcript = getAttributeField(as.character(ret$group), identifier, sep)
        ret = split(ret, transcript)
    }
    on.exit(suppressWarnings(close(con)))
    return(ret)
}