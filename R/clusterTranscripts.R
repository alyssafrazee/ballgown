#' group a gene's assembled transcripts into clusters
#'
#' @param gene name of gene whose transcripts will be clustered.  When using 
#'   Cufflinks output, usually of the form \code{"XLOC_######"}
#' @param gown ballgown object containing experimental data
#' @param k number of clusters to use
#' @param method clustering method to use.  Must be one of \code{"hclust"}, for
#'   hierarchical clustering, or \code{"kmeans"}, for k-means clustering. 
#' 
#' @return list with elements \code{clusters} and \code{pctvar}. 
#'    \code{clusters} contains columns "cluster" and "t_id", and denotes which
#'    transcripts belong to which clusters.  
#' \code{pctvar} is only non-NULL when using k-means clustering and is the
#'    percentage of variation explained by these clusters, defined as the ratio
#'    of the between-cluster sum of squares to the total sum of squares.
#' 
#' @seealso \code{\link{hclust}}, \code{\link{kmeans}}, 
#'   \code{\link{plotLatentTranscripts}} for visualizing the transcript clusters
#' 
#' @author Alyssa Frazee
#' 
#' @export
#' 
#' @examples 
#' data(bg)
#' clusterTranscripts('XLOC_000454', bg, k=2, method='kmeans')
#' # transcripts 1294 and 1301 cluster together, 91% variation explained.

clusterTranscripts = function(gene, gown, k=NULL, method=c('hclust', 'kmeans')){
    method = match.arg(method)

    txnames = indexes(gown)$t2g$t_id[indexes(gown)$t2g$g_id == gene]
    if(!is.null(k) & length(txnames) <= k){
        stop('k must be strictly less than the number of transcripts in gene')
    }else if(is.null(k)){
        k = ceiling(sqrt(length(txnames)/2))
    }

    inds = which(names(structure(gown)$trans) %in% txnames)
    tx = structure(gown)$trans[inds]
    chr = texpr(gown, 'all')$chr[inds[1]]

    xg = expand.grid(c(1:length(tx)), c(1:length(tx)))
    matInds = as.vector(t(xg)) 
    #^ vector of transcript pairs to calculate overlaps for
    olGroup = rep(1:(length(matInds)/2), each=2)
    olGroupSplit = rep(olGroup, times=elementNROWS(tx[matInds]))
    overlapping = split(unlist(tx[matInds]), olGroupSplit)
    coverages = coverage(ranges(overlapping))
    runvals = runValue(coverages)
    runlengths = runLength(coverages)
    pct = sum(runlengths[runvals==2]) / sum(runlengths[runvals==1 | runvals==2])
    distMat = 1-matrix(pct, nrow=length(tx))


    if(method=="hclust"){
        groups = cutree(hclust(as.dist(distMat)), k=k)
        pctvar = NULL
    } 
    
    if(method=="kmeans"){
        km = kmeans(distMat, centers=k)
        groups = as.numeric(km$cluster)
        pctvar = km$betweenss/km$totss
    }

    return(list(clusters=data.frame(cluster=groups, t_id=txnames), 
        pctvar=pctvar))
}
