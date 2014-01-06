#' group a gene's assembled transcripts into clusters
#'
#' @param gene name of gene whose transcripts will be clustered.  When using Cufflinks output, usually of the form \code{"XLOC_######"}
#' @param gown ballgown object containing experimental data
#' @param k number of clusters to use
#' @param method clustering method to use.  Must be one of \code{"hclust"}, for hierarchical clustering, or \code{"kmeans"}, for k-means clustering. 
#' @return list with elements \code{clusters} and \code{pctvar}.  \code{clusters} contains columns "cluster", "tname", and "tid", and denotes which transcripts belong to which clusters.  \code{pctvar} is only non-NULL when using k-means clustering and is the percentage of variation explained by these clusters, defined as the ratio of the between-cluster sum of squares to the total sum of squares.
#' @seealso \code{\link{hclust}}, \code{\link{kmeans}}, \code{\link{\plotLatentTranscripts}} for visualizing the transcript clusters
#' @export
clusterTranscripts = function(gene, gown, k=NULL, method=c("hclust", "kmeans")){
    method = match.arg(method)

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
