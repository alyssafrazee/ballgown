#' cluster a gene's transcripts and calculate cluster-level expression
#'
#' @param gown ballgown object containing experimental data
#' @param dattype which transcript-level expression measurement to use (\code{'cov'}, average 
#' per-base coverage, or \code{'FPKM'})
#' @param method which clustering method to use (\code{'hclust'}, hierarchical clustering, or 
#' \code{'kmeans'}, k-means clustering)
#' 
#' @return data frame with one row per transcript cluster and one column per sample, where entries 
#' are summed expression measurements for all the transcripts in the appropriate cluster
#' 
#' @details Transcript clustering methods are in development, so use this function with caution: for
#' example, it's not clear that the appropriate cluster-level expression measurement is the sum.  
#' 
#' Also, this function runs clustering and collapsing on the entire ballgown object and could be 
#' very slow, so you may want to check out the \code{\link{subset}} method for ballgown objects and 
#' run this function on small chunks of genes.
#' 
#' @seealso \code{\link{hclust}}, \code{\link{kmeans}}, \code{\link{clusterTranscripts}} for 
#' gene-level transcript clustering, \code{\link{plotLatentTranscripts}} for visualizing transcript 
#' clusters
#' 
#' @author Alyssa Frazee
#' 
#' @export
collapseTranscripts = function(gown, dattype=c('cov','FPKM'), method=c('hclust', 'kmeans')){
    dattype = match.arg(dattype)
    method = match.arg(method)
  
    # number of transcripts for each gene:
    txnums = table(expr(gown)$trans$gene_id)
    genes.uniq = names(txnums)
  
    # number of clusters for each transcript:
    k = sapply(genes.uniq, function(g) ifelse(txnums[[g]]<4, txnums[[g]], 
        ceiling(sqrt(txnums[[g]]/2)))) 
    #^^number of clusters for each gene - only cluster if >3 transcripts
  
    # set up the tx-by-sample table:
    tab = matrix(NA, nrow=sum(k), ncol=nrow(indexes(gown)$pData))
  
    # can test either "cov" or "FPKM" (need to choose one)
    coltypes = sapply(names(expr(gown)$trans), gettype)
    columns = which(coltypes==dattype)
  
    # set up vector to hold names of your clusters:
    tclustnames = NULL
    
    # for each gene, cluster transcripts & sum counts to fill table:
    for(g in genes.uniq){
        ind = which(genes.uniq == g)
        if(ind==1) tabinds = c(1:(cumsum(k)[1]))
        if(ind!=1) tabinds = c((cumsum(k)[ind-1]+1):(cumsum(k)[ind]))
        txdata = expr(gown)$trans[expr(gown)$trans$gene_id==g,]
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
            tids = as.numeric(sapply(as.character(clusters[[i]]), function(x){
                substr(x, 3, nchar(x))
            })) 
            trows = which(txdata$t_id %in% tids)
            sumdat[i,] = apply(txdata[trows, columns], 2, sum)
        }
        tab[tabinds,] <- sumdat
        tclustnames[tabinds] <- as.character(sapply(clusters, function(x) paste(x, collapse="-")))
        }
    }
  
    tab = as.data.frame(tab)
    names(tab) = as.character(sapply(names(expr(gown)$trans)[coltypes==dattype], getsamp))
    rownames(tab) = tclustnames
    return(tab)
}  
