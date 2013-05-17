clusterTranscripts = function(gene, gown, k=NULL, method="hclust"){
  txnames = indexes(gown)$t2g$t_id[indexes(gown)$t2g$g_id == gene]
  strucnames = as.numeric(substr(names(structure(gown)$trans),3,nchar(names(structure(gown)$trans))))
  inds = which(strucnames %in% txnames)
  tx = structure(gown)$trans[inds]
  chr = data(gown)$trans$chr[inds[1]] #add error check later, in case the tx's are on different chromosomes.  ugh.
  
  covind = unique(as.numeric(lapply(tx, function(x) runValue(seqnames(x)))))
  covs = lapply(tx, function(x) coverage(x)[[covind]])
  
  startpos = unlist(lapply(covs, function(x) runLength(x)[1]+1))
  endpos = unlist(lapply(covs, function(x) sum(runLength(x))))
  minbp = min(startpos)
  maxbp = max(endpos)
  
  compactrles = lapply(covs, function(x) Rle(values = runValue(x)[-1], lengths = runLength(x)[-1]))
  expanded = lapply(compactrles, IRanges::as.vector)
  tnames = names(expanded)
  expanded = lapply(1:length(expanded), function(i){
    out = expanded[[i]]
    if(startpos[i] > minbp){
      out = c(rep(0, startpos[i]-minbp), out)
    }
    if(endpos[i] < maxbp){
      out = c(out, rep(0, maxbp-endpos[i]))
    }
    return(out)
  })
  d = dist(matrix(unlist(expanded), nrow=length(expanded), byrow=TRUE))
  if(is.null(k)) k = ceiling(sqrt(length(tnames)/2))
  
  if(method=="hclust"){
    h = hclust(d)
    groups = cutree(h, k=k)
    pctvar = NULL
  }
  
  if(method=="kmeans"){
    km = kmeans(as.matrix(d), centers=k)
    groups = as.numeric(km$cluster)
    pctvar = km$betweenss/km$totss
  }
  
  return(list(clusters = data.frame(cluster=groups, tname=tnames, tid=txnames), pctvar = pctvar))
}