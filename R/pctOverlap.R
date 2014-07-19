#' calculate percent overlap between two GRanges objects
#'
#' @param tx1 GRanges object
#' @param tx2 GRanges object
#' 
#' @return percent overlap between \code{tx1} and \code{tx2}, as defined by the 
#'   ratio of the intersection of \code{tx1} and \code{tx2} to the union of 
#'   \code{tx1} and \code{tx2}. 
#' 
#' @details In the ballgown context, \code{tx1} and \code{tx2} are two
#'   transcripts, each represented by GRanges objects whose ranges represent the
#'   exons comprising the transcripts.  The percent overlap is the number of 
#'   nucleotides falling within both transcripts divided by the number of 
#'   nucleotides falling within either transcript.  Useful as a measure of 
#'   transcript closeness (as it is essentially Jaccard distance).

#' @author Alyssa Frazee

#' @export
#' 
#' @examples
#' data(bg)
#' gtfPath = system.file('extdata', 'annot.gtf.gz', package='ballgown')
#' annot_grl = gffReadGR(gtfPath, splitByTranscript=TRUE)
#' pctOverlap(structure(bg)$trans[[2]], annot_grl[[369]]) #79.9%
pctOverlap = function(tx1, tx2){
    stopifnot(class(tx1) == 'GRanges' & class(tx2) == 'GRanges')
    ch1 = as.character(runValue(seqnames(tx1)))
    ch2 = as.character(runValue(seqnames(tx2)))
    if(ch1 != ch2){
        return(0)
    }
    tmp1 = reduce(tx1)
    tmp2 = reduce(tx2)
    if(!identical(ranges(tmp1), ranges(tx1))){
        warning('tx1 contained overlapping ranges and was reduced.')
        tx1 = tmp1
    }
    if(!identical(ranges(tmp2), ranges(tmp2))){
        warning('tx2 contained overlapping ranges and was reduced.')
        tx2 = tmp2
    }
    mcols(tx1) = NULL
    mcols(tx2) = NULL # take away metadata
    ntcov = coverage(c(tx1, tx2))
    ind = which(names(ntcov)==ch1)
    covrle = ntcov[[ind]]
    covrle_val = runValue(covrle)
    covrle_len = runLength(covrle)
    return(sum(covrle_len[covrle_val==2]) / sum(covrle_len[covrle_val==1 | 
        covrle_val == 2]))
}
