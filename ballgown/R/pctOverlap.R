# function to calculate percent overlap between 2 GRanges objects

pctOverlap = function(tx1, tx2){
    ch1 = as.character(runValue(seqnames(tx1)))
    ch2 = as.character(runValue(seqnames(tx2)))
	stopifnot(ch1 == ch2)
	ntcov = coverage(c(tx1, tx2))
	ind = which(names(ntcov)==thechr)
	covrle = ntcov[[ind]]
	covrle_val = runValue(covrle)
	covrle_len = runLength(covrle)
	return(sum(covrle_len[covrle_val==2]) / sum(covrle_len[covrle_val==1 | covrle_val == 2]))
}
