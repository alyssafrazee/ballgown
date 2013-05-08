plotLatentTranscripts = function(gene, gown, km = NULL, k = NULL, userefseq = TRUE, gtf = NULL, genome="hg19", position="txStart", returnkm = TRUE, ...){
	
	### ... are extra arguments for getRefSeq
	
	# work with reasonably small gown object
	cond = paste0("gene_id == \"", gene, "\"")
	gown = subset(gown, cond, global=FALSE)

	require(RColorBrewer)
	
	# plot setup:
	ma = IRanges::as.data.frame(structure(gown)$trans)
	thetranscripts = indexes(gown)$t2g$t_id[indexes(gown)$t2g$g_id==gene]
	thetranscripts = paste0("tx", thetranscripts)
	gtrans = subset(ma, element %in% thetranscripts)
#	gtrans$tid = as.numeric(sapply(gtrans$element, function(x) as.numeric(substr(x,3,nchar(x))))) #why?
	
	if(userefseq){
		refinfo = getRefSeq(gene, gown, gtf, genome = genome, position = position)
#		refinfo$gtf$tid = refinfo$gtf$element  #???
		gtrans = rbind(gtrans, IRanges::as.data.frame(refinfo$gtf))
	}
	
	xax = seq(min(gtrans$start), max(gtrans$end), by=1)
	
	# these might mess w/ refseq stuff.
	# if(length(unique(gtrans$seqnames)) > 1) stop("Your gene appears to span multiple chromosomes, which is interesting but also kind of annoying, R-wise.  Please choose another gene until additional functionality is added!")
	# if(length(unique(gtrans$strand)) > 1) stop("Your gene appears to contain exons from both strands, which is potentially interesting but also kind of confusing, so please choose another gene until we figure this sucker out.")


	if(is.null(km)){
		dmat = 100*(1-transcriptOverlaps(gene, gown, userefseq = userefseq, gtf))
		# option 1: do k-means clustering with automatic choice of k
		if(is.null(k)){
			for(i in 1:(nrow(dmat)-1)){
				km = kmeans(dmat, centers=i)
				pctvar = km$betweenss/km$totss
				if(pctvar>=0.9) break
			}
			if(i==(nrow(dmat)-1)){
				# this means nothing explained 90% of the variation
				plotTranscripts(gene, samp = paste0("cov.",as.character(indexes(gown)$pData$dirname[1])), gown, legend = TRUE, colorby="transcript")
				warning("best choice of k is approximately the # of assembled transcripts. we recommend either not clustering these at all, or using the # of RefSeq transcripts for k")
			}
		}
		# option 2: do k-means clustering with fixed k
		if(!is.null(k)){
			km = kmeans(dmat, centers = k)
		}
	}
	
	# PLOT:
	plot(xax, rep(0,length(xax)), ylim=c(0,nrow(dmat)+1), type="n", xlab="genomic position", yaxt = "n", ylab="")

	par(mar=c(5,2,4,2))
	
	
	cols = brewer.pal(length(unique(km$cluster)), "Dark2")
	for(tx in names(sort(km$cluster))){
		txind = which(names(sort(km$cluster))==tx)
		gtsub = gtrans[gtrans$element==tx,]
		gtsub = gtsub[order(gtsub$start),]
		for(exind in 1:dim(gtsub)[1]){
			mycolor = ifelse(substr(tx,1,2)=="tx", cols[sort(km$cluster)[txind]], "gray70")
			bcolor = ifelse(substr(tx,1,2)=="tx", "black", cols[sort(km$cluster)[txind]])
			polygon(x=c(gtsub$start[exind], gtsub$start[exind], gtsub$end[exind], gtsub$end[exind]), y=c(txind-0.4, txind+0.4, txind+0.4, txind-0.4), col=mycolor, border=bcolor)
			if(exind != dim(gtsub)[1]){
				lines(c(gtsub$end[exind], gtsub$start[exind+1]), c(txind, txind), lty=2, col="gray60")
			}
		}
	}
	
	k = length(km$size)
	title(paste0(gene,": transcripts clustered with k-means, k=",k))
	if(returnkm){
		return(km)
	}
}
