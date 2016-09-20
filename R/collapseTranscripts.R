#' cluster a gene's transcripts and calculate cluster-level expression
#'
#' @param gene which gene's transcripts should be clustered
#' @param gown ballgown object
#' @param meas which transcript-level expression measurement to use 
#'   (\code{'cov'}, average per-base coverage, or \code{'FPKM'})
#' @param method which clustering method to use: \code{'hclust'} (hierarchical
#'   clustering) or \code{'kmeans'} (k-means clustering).
#' @param k how many clusters to use. 
#' 
#' @return list with two elements:
#' \itemize{
#'   \item \code{tab}, a cluster-by-sample table of expression measurements
#'     (\code{meas}, either cov or FPKM), where the expression measurement for 
#'     each cluster is the mean (for \code{'cov'}) or aggregate (for 
#'     \code{'FPKM'}, as in \code{\link{gexpr}}) expression measurement for all
#'     the transcripts in that cluster. This table can be used as the 
#'     \code{gowntable} argument to \code{\link{stattest}}, if differential
#'     expression results for transcript *clusters* are desired. 
#'   \item \code{cl} output from \code{\link{clusterTranscripts}} that was run 
#'     to produce \code{tab}, for reference. Cluster IDs in the \code{cluster} 
#'     component correspond to row names of \code{tab}
#' }
#' 
#' @seealso \code{\link{hclust}}, \code{\link{kmeans}}, 
#'    \code{\link{clusterTranscripts}}, \code{\link{plotLatentTranscripts}}
#' 
#' @author Alyssa Frazee
#' 
#' @export
#' @examples 
#' data(bg)
#' collapseTranscripts(bg, gene='XLOC_000454', meas='FPKM', method='kmeans')
#' 
collapseTranscripts = function(gene, gown, meas='FPKM', 
    method=c('hclust', 'kmeans'), k=NULL){

    if(!gown@RSEM){
        meas = match.arg(meas, c('cov', 'FPKM'))
    }else{
        meas = match.arg(meas, c('FPKM', 'TPM'))
    }
    method = match.arg(method)
    if(is.null(k)){
        k = 'thumb'
    } else if(is.character(k)){
        k = match.arg('thumb', 'var90')
        if(k == 'var90' & method != 'kmeans'){
            stop(.makepretty('choosing k explaining 90% of variability is only
                available with kmeans clustering.'))
        }
    } else if(as.integer(k) != k){
        k = floor(k)
        warning('k must be an integer. Rounding down.')
    } else if(!is.numeric(k)){
        stop('k must be "var90", "thumb", or an integer (number of clusters)')
    }

    gown = subset(gown, paste0("gene_id == '", gene, '\''))
    ntranscripts = length(structure(gown)$trans)

    # cluster:
    if(k == 'thumb'){
        k = ceiling(sqrt(ntranscripts/2))
        cl = clusterTranscripts(gene=gene, gown=gown, method=method, k=k)
        message(paste('using k =', k, 'clusters'))
    } else if(k == "var90"){
        for(i in 1:(ntranscripts-1)){
            cl = clusterTranscripts(gene=gene, gown=gown, method="kmeans", k=i)
            if(cl$pctvar>=0.9) break
        }
        if(i==ntranscripts-1 & cl$pctvar < 0.9){
            # either k=n-1 or nothing explained 90% of variation:
            stop(.makepretty("k = n-1 did not explain 90% of variation. Try no
                clustering, or specifying k."))
        }
    } else {
        cl = clusterTranscripts(gene=gene, gown=gown, method=method, k=k)
    }

    clusterList = split(cl$clusters$t_id, cl$clusters$cluster)
    collapsed = lapply(clusterList, function(x){
        if(length(x) > 1){
            if(meas == 'cov'){
                colMeans(texpr(gown)[texpr(gown,'all')$t_id %in% x,])    
            } else {
                tstruct = structure(gown)$trans[names(structure(gown)$trans) 
                    %in% x]
                tlengths = sapply(width(tstruct), sum)
                clength = sum(width(reduce(unlist(tstruct))))
                totfrags = colSums(tlengths * texpr(gown)[texpr(gown,'all')$t_id
                    %in% x,])
                totfrags / clength
            }
        }else{
            texpr(gown)[texpr(gown,'all')$t_id == x,]
        }
    })

    tab = t(as.data.frame(collapsed))
    rownames(tab) = names(collapsed)
    return(list(tab=tab, cl=cl))
}
