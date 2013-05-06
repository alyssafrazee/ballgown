getRefSeq = function(gene, gown, gtf, genome = "hg19", position="txStart"){
	cond = paste0("gene_id ==\"", gene,"\"")
	gown = subset(gown, cond, global = FALSE)
	chrom = unique(data(gown)$trans$chr)
	midp = (min(start(ranges(structure(gown)$exon))) + max(end(ranges(structure(gown)$exon))))/2
	if(substr(chrom,1,1)!="c") chrom = paste0("chr",chrom) ### THIS MAY NEED SOME EDITING
	require(ACME)
	gname = as.character(findClosestGene(chrom, pos = midp, genome = genome, position = position )$name)
	gtf.tmp = subset(gtf, (geneid %in% gname) & seqname==chrom )
	isdup = grepl("dup", gtf.tmp$txid)
	gtf.df = gtf.tmp[!isdup,]
	txlist = split(gtf.df, gtf.df$txid)
	require(GenomicRanges)
	txlist = lapply(txlist, function(x) GRanges(seqnames = Rle(x$seqname), ranges=IRanges(start=x$start, end=x$end), strand = x$strand, id = rep(0,length(x$strand)) , transcripts=rep("NA", length(x$strand)) ) )
	ret = list(gene=gname, gtf = GRangesList(txlist))
	return(ret)
	}

