# function to calculate percent overlap between 2 GRanges objects

pctOverlap = function(tx1, tx2){
	thechr = runValue(c(seqnames(tx1), seqnames(tx2)))
	stopifnot(length(thechr) == 1)
	ntcov = coverage(c(tx1, tx2))
	ind = which(names(ntcov)==thechr)
	covrle = ntcov[[ind]]
	covrle_val = runValue(covrle)
	covrle_len = runLength(covrle)
	return(sum(covrle_len[covrle_val==2]) / sum(covrle_len[covrle_val==1 | covrle_val == 2]))
}
