transcriptOverlaps = function(gene, gown, userefseq = TRUE, gtf = NULL, genome = "hg19", position = "txStart"){
	
	require(GenomicRanges)
	
	txnames = indexes(gown)$t2g$t_id[indexes(gown)$t2g$g_id == gene]
	strucnames = as.numeric(substr(names(structure(gown)$trans),3,nchar(names(structure(gown)$trans))))
	inds = which(strucnames %in% txnames)
	tx = structure(gown)$trans[inds]
	
	if(userefseq){
		reftx = getRefSeq(gene, gown, gtf, genome, position)$gtf
		txchrs = seqlevels(tx)
		refseqchrs = seqlevels(reftx)
		if( any(!(refseqchrs %in% txchrs)) ){
			# I will assume for now that this means our transcripts are labeled "1", "2", ... and RefSeq is labeled "chr1", "chr2", etc. this can be fixed later.
			seqlevels(tx) = as.character(sapply(seqlevels(tx), function(x) paste0("chr",x)))
		}
		tx = append(tx, reftx)
	}
	
	overlapMat = matrix(NA, nrow=length(tx), ncol=length(tx))
	rownames(overlapMat) = colnames(overlapMat) = names(tx)
	diag(overlapMat) = 1
	for(ii in 1:length(tx)){
		for(jj in 1:length(tx)){
			if(ii==jj) next
			bpsi = bpset(tx[[ii]])
			bpsj = bpset(tx[[jj]])
			overlapMat[ii,jj] = length(intersect(bpsi, bpsj))/length(bpsi)
		}
	}
	
	return(overlapMat)
}  #number in row i, column j answers the question "what percentage of transcript i is overlapped by transcript j?"
