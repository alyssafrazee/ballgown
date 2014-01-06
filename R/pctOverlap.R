#' calculate percent overlap between two GRanges objects
#'
#' @param tx1 GRanges object
#' @param tx2 GRanges object
#' @return percent overlap between \code{tx1} and \code{tx2}, as defined by the ratio of the intersection of \code{tx1} and \code{tx2} to the union of \code{tx1} and \code{tx2}.
#' @details In the ballgown context, \code{tx1} and \code{tx2} are two transcripts, each represented by GRanges objects whose ranges represent the exons comprising the transcripts.  The percent overlap is the number of nucleotides falling within both transcripts divided by the number of nucleotides falling within either transcript.  Useful as a measure of transcript closeness. 
#' @author Alyssa Frazee
#' @export
pctOverlap = function(tx1, tx2){
    ch1 = as.character(runValue(seqnames(tx1)))
    ch2 = as.character(runValue(seqnames(tx2)))
	stopifnot(ch1 == ch2)
	ntcov = coverage(c(tx1, tx2))
	ind = which(names(ntcov)==ch1)
	covrle = ntcov[[ind]]
	covrle_val = runValue(covrle)
	covrle_len = runLength(covrle)
	return(sum(covrle_len[covrle_val==2]) / sum(covrle_len[covrle_val==1 | covrle_val == 2]))
}
