difftest = function(gown, feature = "transcript", method = "DESeq", dattype = "cov", collapse = TRUE, ...){
	if(feature=="transcript"){
		if(dattype!="cov" & dattype!="FPKM") stop("transcripts only have cov and FPKM measurements")
		coltypes = as.character(sapply(names(data(gown)$trans), gettype))
		inds = which(coltypes==dattype)
		
		if(collapse){
			# merge the transcripts using k-means, reduce, and sum function
			# this will be super slow at first, since we have to do merging for all the genes...this overlap stuff needs to be way faster.
			for(g in unique(data(gown)$trans$gene_id)){
				dmat = 100*(1-transcriptOverlaps(gene, gown, userefseq = userefseq, gtf))
				## choose k:
				for(i in 1:(nrow(dmat)-1)){
					km = kmeans(dmat, centers=i)
					pctvar = km$betweenss/km$totss
					if(pctvar>=0.9) break
				}
				if(i==(nrow(dmat)-1)){
					km = kmeans(dmat, centers = nrow(dmat)-1)
					warning(paste0("gene ", g,": number of clusters is approx. number of transcripts")
				}
				
				

		}
		
		
		tab = data(gown)$trans[,inds]
	}
	if(feature=="exon"){
		exontypes = unique(as.character(sapply(names(data(gown)$exon)[-c(1:5)], gettype)))
		if(!(dattype %in% exontypes)) stop(paste0("exons only have the following measurements: ", paste(exontypes, collapse=", ")))
		coltypes = as.character(sapply(names(data(gown)$exon), gettype))
		inds = which(coltypes==dattype)
		tab = data(gown)$exon[,inds]
	}
	if(feature=="junction" | feature=="intron"){
		introntypes = c("rcount", "ucount", "mrcount")
		if(!(dattype %in% introntypes)) stop(paste0("introns/junctions only have the following measurements: ", paste(introntypes, collapse=", ")))
		coltypes = as.character(sapply(names(data(gown)$intron), gettype))
		inds = which(coltypes==dattype)
		tab = data(gown)$intron[,inds]
	}
	if(!(feature %in% c("transcript","exon","junction","intron"))) stop(paste0("unknown feature: ",feature,". Please choose one of transcript, exon, intron, or junction.  (\"intron\" and \"junction\" are the same test)."))
	
	
	
}