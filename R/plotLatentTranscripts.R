#' cluster assembled transcripts and plot the results
#' 
#' This is an experimental, first-pass function that clusters assembled 
#' transcripts based on their overlap percentage, then plots and colors the 
#' transcript clusters.
#'
#' @param gene string, name of gene whose transcripts should be clustered 
#'   (e.g., "XLOC_000001")
#' @param gown object of class \code{ballgown} being used for analysis
#' @param method clustering method to use.  Currently can choose from 
#'   hierarchical clustering (\code{hclust}) or K-means (\code{kmeans}).  More 
#'   methods are in development.
#' @param k number of transcripts clusters to use.  By default, \code{k} is 
#'   \code{NULL} and thus is chosen using a rule of thumb, but providing
#'   \code{k} overrides those rules of thumb.
#' @param choosek if \code{k} is not provided, how should the number of clusters
#'   be chosen?  Must be one of "var90" (choose a \code{k} that explains 90 
#'   percent of the observed variation) or "thumb" (\code{k} is set to be 
#'   approximately \code{sqrt(n)}, where n is the total number of transcripts 
#'   for \code{gene})
#' @param returncluster if TRUE (as it is by default), return the results of the
#'   call to \code{clusterTrancsripts} so the data is available for later use.  
#'   Nothing is returned if FALSE.
#' @param labelTranscripts if TRUE (as it is by default), print transcript IDs 
#'   on the y-axis
#' @param ... other arguments to pass to plotTranscripts
#' 
#' @return if \code{returncluster} is TRUE, the transcript clusters are returned
#'   as described in \code{\link{clusterTranscripts}}. A plot of the transcript 
#'   clusters is also produced, in the style of \code{\link{plotTranscripts}}.
#' 
#' @seealso \code{\link{clusterTranscripts}}, \code{\link{plotTranscripts}}
#' 
#' @author Alyssa Frazee
#' 
#' @export
#' @examples \donttest{
#' data(bg)
#' plotLatentTranscripts('XLOC_000454', bg, method='kmeans', k=2)
#' }
plotLatentTranscripts = function(gene, gown, method=c("hclust", "kmeans"), 
    k=NULL, choosek=c("var90", "thumb"), returncluster=TRUE, 
    labelTranscripts=TRUE, ...){

    ## check validity:
    method = match.arg(method)
    if(is.null(k)){
        choosek = match.arg(choosek)
        if(choosek == 'var90' & method != 'kmeans'){
            stop("need to use method=\"kmeans\" when choosek is \"var90\"")
        }
    }

    ## do the clustering
    if(!is.null(k)){
        cl = clusterTranscripts(gene=gene, gown=gown, method=method, k=k)
    } else if(choosek=="thumb"){
        ntranscripts = sum(texpr(gown, 'all')$gene_id == gene)
        k = ceiling(sqrt(ntranscripts/2))
        cl = clusterTranscripts(gene=gene, gown=gown, method=method, k=k)
    } else if(choosek == "var90"){
        ntranscripts = sum(texpr(gown, 'all')$gene_id == gene)
        for(i in 1:(ntranscripts-1)){
            cl = clusterTranscripts(gene=gene, gown=gown, method="kmeans", k=i)
            if(cl$pctvar>=0.9) break
        }
        if(i==ntranscripts-1 & cl$pctvar < 0.9){
            # either k=n-1 or nothing explained 90% of variation:
            plotTranscripts(gene=gene, gown=gown, colorby="none")
            warning(.makepretty("k = n-1 did not explain 90% of variation. Try
                no clustering, or specifying k."))
            return(NULL) #this function is now done.    
        }
    }
    
    ## plot:
    numcolors = length(unique(cl$clusters$cluster))
    cols = suppressWarnings(brewer.pal(numcolors, 'Dark2')) 
    #^have 2-color graph w/o warning
    sorted_clusters = cl$clusters[order(cl$clusters$cluster),]
    pTitle = paste0(gene, ": transcripts clustered with ", method, 
        ", k=", numcolors)
    plotTranscripts(gene, gown, legend=FALSE, samples=sampleNames(gown)[1], 
        main=pTitle, customCol=cols[sorted_clusters$cluster], 
        customOrder=sorted_clusters$t_id, labelTranscripts=labelTranscripts, 
        ...)

    if(returncluster) return(cl) 
}
