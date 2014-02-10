#' read in gtf file as GRanges object
#' 
#' (very) light wrapper for rtracklayer::import
#'  
#' @param gtf name of GTF/GFF file on disk
#' @return object of class \code{GRanges} representing the genomic features in \code{gtf}
#' @seealso \code{\link{gffRead}} for reading in a GTF file as a data frame rather than GRanges.
#' @author Alyssa Frazee
#' @export
gffReadGR = function(gtf){
    require(rtracklayer)
    import(file(gtf), format='GFF')
}