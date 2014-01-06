collapseTranscripts = function(gown, dattype=c('cov','FPKM'), method=c('hclust', 'kmeans')){
  dattype = match.arg(dattype)
  method = match.arg(method)
  
  # number of transcripts for each gene:
  txnums = table(data(gown)$trans$gene_id)
  genes.uniq = names(txnums)
  
  # number of clusters for each transcript:
  k = sapply(genes.uniq, function(g) ifelse(txnums[[g]]<4, txnums[[g]], ceiling(sqrt(txnums[[g]]/2)))) #number of clusters for each gene - only cluster if >3 transcripts
  
  # set up the tx-by-sample table:
  tab = matrix(NA, nrow=sum(k), ncol=nrow(indexes(gown)$pData))
  
  # can test either "cov" or "FPKM" (need to choose one)
  coltypes = sapply(names(data(gown)$trans), gettype)
  columns = which(coltypes==dattype)
  
  # set up vector to hold names of your clusters:
  tclustnames = NULL
    
  # for each gene, cluster transcripts & sum counts to fill table:
  for(g in genes.uniq){
    ind = which(genes.uniq == g)
    if(ind==1) tabinds = c(1:(cumsum(k)[1]))
    if(ind!=1) tabinds = c((cumsum(k)[ind-1]+1):(cumsum(k)[ind]))
    txdata = data(gown)$trans[data(gown)$trans$gene_id==g,]
    if(k[ind]==txnums[[ind]]){
      tab[tabinds,]  <- as.matrix(txdata[, columns])
      tclustnames[tabinds] <- paste0("tx", txdata$t_id)
    }
    if(k[ind] != txnums[[ind]]){
      # do clustering:
      cl = clusterTranscripts(gene=g, gown=gown, k=k[ind], method=method)
      clusters = split(cl$clusters$tname, cl$clusters$cluster)
      sumdat = matrix(NA, nrow=k[ind], ncol=length(columns))
      for(i in 1:(k[ind])){
        tids = as.numeric(sapply(as.character(clusters[[i]]), function(x) substr(x, 3, nchar(x))))
        trows = which(txdata$t_id %in% tids)
        sumdat[i,] = apply(txdata[trows, columns], 2, sum)
      }
      tab[tabinds,] <- sumdat
      tclustnames[tabinds] <- as.character(sapply(clusters, function(x) paste(x, collapse="-")))
    }
  }
  
  tab = as.data.frame(tab)
  names(tab) = as.character(sapply(names(data(simgown)$trans)[coltypes==dattype], getsamp))
  rownames(tab) = tclustnames
  return(tab)
}  
