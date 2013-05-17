plotLatentTranscripts = function(gene, gown, k = NULL, choosek = c("var90", "thumb"), returncluster = TRUE, method=c("hclust", "kmeans")){
  
  require(RColorBrewer)
  
  ## check validity:
  if(length(choosek)==2){
    choosek = "preset"
  }
  if(choosek=="var90" & method!="kmeans"){
    stop("need to use method=\"kmeans\" when choosek is \"var90\"")
  }
  if(is.null(k) & choosek=="preset"){
    stop("must specify either k or a method of choosing k (var90 or thumb)")
  }
  if(length(method)==2){
    stop("please specify a method (hclust or kmeans)")
  }
  
	## work with reasonably small gown object
	cond = paste0("gene_id == \"", gene, "\"")
	gown = subset(gown, cond, global=FALSE)
	
	## plot setup:
	ma = IRanges::as.data.frame(structure(gown)$trans)
	thetranscripts = indexes(gown)$t2g$t_id[indexes(gown)$t2g$g_id==gene]
	thetranscripts = paste0("tx", thetranscripts)
	gtrans = subset(ma, element %in% thetranscripts)
	xax = seq(min(gtrans$start), max(gtrans$end), by=1)
	
  ## check validity again:
	if(length(unique(gtrans$seqnames)) > 1) stop("Your gene appears to span multiple chromosomes, which is interesting but also kind of annoying, R-wise.  Please choose another gene until additional functionality is added!")
	if(length(unique(gtrans$strand)) > 1) stop("Your gene appears to contain exons from both strands, which is potentially interesting but also kind of confusing, so please choose another gene until we figure this sucker out.")
  
  ## do the clustering
  if(!is.null(k)){
    cl = clusterTranscripts(gene=gene, gown=gown, method=method, k=k)
  }
  if(choosek=="thumb"){
    k = ceiling(sqrt(length(thetranscripts)/2))
    cl = clusterTranscripts(gene, gown=gown, method=method, k=k)
  }
  if(choosek == "var90"){
    for(i in 1:(length(thetranscripts)-1)){
      cl = clusterTranscripts(gene=gene, gown=gown, method="kmeans", k=i)
      if(cl$pctvar>=0.9) break
    }
    if(i==length(thetranscripts)-1){
      # nothing explained 90% of variation:
      plotTranscripts(gene=gene, gown=gown, samp=NULL, colorby="none")
      warning("k = n-1 did not explain 90% of variation. we recommend not clustering this gene or pre-specifying a reasonable k for this gene.")
      return(NULL) #this function is now done.
    }
  }
  
	
	## PLOT:
  par(mar=c(5,2,4,2))
  plot(xax, rep(0,length(xax)), ylim=c(0,length(thetranscripts)+1), type="n", xlab="genomic position", yaxt = "n", ylab="")
  
  cols = suppressWarnings(brewer.pal(length(unique(cl$clusters$cluster)), "Dark2"))
  # suppress warnings because we want to be able to have a 2-color graph w/o a warning
  
  clusters.sorted = cl$clusters[order(cl$clusters$cluster),]
  tx.new = clusters.sorted$tname
  cid = clusters.sorted$cluster
  
	for(tx in tx.new){
		txind = which(tx.new==tx)
		gtsub = gtrans[gtrans$element==tx,]
		gtsub = gtsub[order(gtsub$start),]
		for(exind in 1:dim(gtsub)[1]){
			mycolor = cols[cid[txind]]
			polygon(x=c(gtsub$start[exind], gtsub$start[exind], gtsub$end[exind], gtsub$end[exind]), y=c(txind-0.4, txind+0.4, txind+0.4, txind-0.4), col=mycolor)
			if(exind != dim(gtsub)[1]){
				lines(c(gtsub$end[exind], gtsub$start[exind+1]), c(txind, txind), lty=2, col="gray60")
			}
		}
	}
	
  numclust = length(unique(cl$clusters$cluster))
	title(paste0(gene,": transcripts clustered with ",method,", k=",numclust))
	if(returncluster){
		return(cl)
	}
}
