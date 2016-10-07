#' visualize transcript abundance by group
#'
#' @param gene name of gene whose transcripts will be plotted.  When using 
#'   Cufflinks/Tablemaker output, usually of the form \code{"XLOC_######"}
#' @param gown ballgown object containing experimental and phenotype data
#' @param overall if \code{TRUE}, color features by the overall 
#'   (experiment-wide) mean rather than a group-specific mean
#' @param groupvar string representing the name of the variable denoting which 
#'   sample belongs to which group.  Can be \code{"none"} (if you want the 
#'   study-wide mean), or must correspond to the name of a column of 
#'   \code{pData(gown)}.  Usually a categorical variable.
#' @param groupname string representing which group's expression means you want 
#'   to plot.  Can be \code{"none"} (if you want the study-wide mean), 
#'   \code{"all"} (if you want a multipanel plot of each group's mean 
#'   expression), or any of the levels of \code{groupvar}.
#' @param meas type of expression measurement to plot. One of "cov", "FPKM", 
#'   "rcount", "ucount", "mrcount", or "mcov". Not all types are valid for all 
#'   features. (See description of tablemaker output for more information).
#' @param colorby one of \code{"transcript"} or \code{"exon"}, indicating which 
#'   feature's abundances should dictate plot coloring. 
#' @param legend if \code{TRUE} (as it is by default), a color legend is drawn 
#'   on top of the plot indicating the scale for feature abundances.
#' @param labelTranscripts if \code{TRUE}, transcript ids are labeled on the 
#'   left side of the plot. Default \code{FALSE}.
#' 
#' @return produces a plot of the transcript structure for the specified gene in
#'   the current graphics device, colored by study-wide or group-specific mean 
#'   expression level.
#' 
#' @seealso \code{\link{plotTranscripts}} 
#' 
#' @author Alyssa Frazee
#' 
#' @export
#' @examples \donttest{
#' data(bg)
#' plotMeans('XLOC_000454', bg, groupvar='group', meas='FPKM', 
#'   colorby='transcript')
#' }

plotMeans = function(gene, gown, overall=FALSE, groupvar, groupname='all', 
    meas=c('cov', 'FPKM', 'rcount', 'ucount', 'mrcount', 'mcov'),
    colorby=c('transcript', 'exon'), legend=TRUE, labelTranscripts=FALSE){
    
    meas = match.arg(meas)
    colorby = match.arg(colorby)

    if(!(meas %in% gown@meas | gown@meas=='all')){
        stop(paste('gown does not include', meas, 'measurements'))
    }

    if(colorby == "transcript" & !(meas %in% c("cov", "FPKM"))){
        stop("transcripts only have cov and FPKM measurements")
    }
    if(colorby == "exon" & meas=="FPKM"){
        stop("exons do not have FPKM measurements")
    }
    if(class(gown)!="ballgown") stop("gown must be a ballgown object")

    if(!overall & (groupvar=="none"|groupname=="none")){
        stop(.makepretty("to plot means for a specific group, please provide
            both the grouping variable name and the specific group name."))
    }
    if(groupname == "all" & overall){
        warning(.makepretty("Plotting the study-wide mean for each feature. To
            plot each group's mean separately, set overall=FALSE."))
    }
    if(groupname == "all" & groupvar == "none"){
        stop("to plot means for all groups, please provide groupvar (the name
            of the group variable)")
    }
    
    # some setup:
    ma = IRanges::as.data.frame(structure(gown)$trans)
    if(names(ma)[2] != 'group_name'){
        stop('IRanges::as.data.frame has changed. Please report as issue.')
    }
    
    thetranscripts = indexes(gown)$t2g$t_id[indexes(gown)$t2g$g_id==gene]
    gtrans = ma[ma$group_name %in% thetranscripts,]
    gtrans$tid = gtrans$group_name
    xax = seq(min(gtrans$start), max(gtrans$end), by=1)
    numtx = length(unique(thetranscripts))
    if(groupvar != "none"){
        pdatacol = which(names(pData(gown))==groupvar) 
    }
    if(groupname == "all" & !overall){
        numplots = length(unique(pData(gown)[,pdatacol]))
        mf = c(floor(sqrt(numplots)), ceiling(sqrt(numplots)))
        plot_titles = sort(unique(as.factor(pData(gown)[,pdatacol])))
    }else{
        numplots = 1
        mf = c(1,1)
        if(overall){
            plot_titles = 'overall'
        }else{
            plot_titles = ifelse(groupname=="none", "mean across groups", 
                groupname)
        }
    }
    if(overall){
        samples = list(sampleNames(gown))
    }else if(groupname == "all"){
        samples = split(sampleNames(gown), pData(gown)[,pdatacol])
    }else{
        samples = list(sampleNames(gown)[pData(gown)[,pdatacol]==groupname]) 
    }
    
    # plot base:
    westval = ifelse(labelTranscripts, 4, 2)
    par(mar=c(5, westval, 4, 2), mfrow = mf)
    ymax = ifelse(legend, numtx+2, numtx+1)

    # the for-loop only has >1 iteration if groupname == "all"
    for(p in 1:numplots){
        plot(xax, rep(0,length(xax)), ylim=c(0,ymax), 
        type="n", xlab="genomic position", 
        main=paste0(gene,": ",plot_titles[p]), 
        yaxt = "n", ylab="")

        # calculate means
        if(colorby == "transcript"){
            # mean transcript-level expression measurement for this group
            t_id = texpr(gown, 'all')$t_id
            expr_matrix_i = which(t_id %in% gtrans$tid)
            smalldat = texpr(gown, meas)[expr_matrix_i,]
            datacols = which(sampleNames(gown) %in% samples[[p]])
            if(length(expr_matrix_i) > 1) {
                tmeans = rowMeans(smalldat[,datacols])
            } else {
                tmeans = mean(smalldat[datacols])
                names(tmeans) = texpr(gown, 'all')$t_id[expr_matrix_i]  # i has length 1
            }
        }

        if(colorby == "exon"){
            # mean exon-level expression measurement for this group
            e_id = eexpr(gown, 'all')$e_id
            smalldat = as.matrix(eexpr(gown, meas)[e_id %in% gtrans$id,])
            if(class(smalldat) != 'matrix'){
                smalldat = t(as.matrix(smalldat))
            }
            datacols = which(sampleNames(gown) %in% samples[[p]])
            emeans = rowMeans(smalldat[,datacols])
        }
    
        # make color scale:
        maxcol = quantile(smalldat, 0.99)
        colscale = seq(0, maxcol, length.out=200)

        # draw the transcripts
        for(tx in unique(gtrans$tid)){
            if(colorby == "transcript") {
                mycolor = closestColor(tmeans[which(names(tmeans)==tx)], colscale)
            }
            txind = which(unique(gtrans$tid)==tx)
            gtsub = gtrans[gtrans$tid==tx,]
            gtsub = gtsub[order(gtsub$start),]
            for(exind in 1:dim(gtsub)[1]){
                if(colorby == "exon"){ 
                    mycolor = closestColor(
                        emeans[which(names(emeans)==gtsub$id[exind])], colscale)
                }
                stopifnot(length(mycolor) > 0)
                polygon(
                    x=c(gtsub$start[exind], gtsub$start[exind], 
                        gtsub$end[exind], gtsub$end[exind]),
                    y=c(txind-0.4,txind+0.4,txind+0.4,txind-0.4), col=mycolor)
                if(exind!=dim(gtsub)[1]){
                    lines(c(gtsub$end[exind],gtsub$start[exind+1]), 
                        c(txind, txind), lty=2, col="gray60")
                }
            }
        }

        # draw the legend:
        if(legend){
            leglocs = seq(min(xax)+1, max(xax)-1, length=length(colscale)+1)
            for(i in 1:length(colscale)){
                polygon(x=c(leglocs[i], leglocs[i], leglocs[i+1], leglocs[i+1]),
                    y=c(ymax-0.3, ymax, ymax, ymax-0.3), border=NA, 
                    col=rev(heat.colors(length(colscale)))[i])
            }
            text(x=seq(min(xax)+1, max(xax)-1, length=20), y=rep(ymax+0.1, 20), 
                labels=round(colscale,2)[seq(1,length(colscale), length=20)], 
                cex=0.5)
            text(x=median(xax), y=ymax-0.5, labels=paste("mean expression, by",
                colorby), cex=0.5)
        }

        if(labelTranscripts){
            axis(side=2, at=c(1:numtx), labels=unique(gtrans$tid), 
                cex.axis=0.75, las=1)
        }
    }
}