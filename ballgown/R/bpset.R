#' get vector of nucleotides comprising a transcript
#'
#' @param tra a GRanges object representing a transcript (see the transcript "structure" element of a ballgown object)
#' @return vector of nucleotide positions included in the transcript
#' @seealso \code{\link{transcriptOverlaps}} which is this function's wrapper
#' @export

bpset = function(tra){
	unique(unlist(sapply(1:length(ranges(tra)), function(i) c(start(tra[i,]):end(tra[i,])))))
}