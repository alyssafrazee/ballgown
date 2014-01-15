#' visualize transcript abundance by group
#'
#' @param gene name of gene whose transcripts will be plotted.  When using Cufflinks output, usually of the form \code{"XLOC_######"}
#' @param gown ballgown object containing experimental and phenotype data
#' @param overall if \code{TRUE}, color features by the overall (experiment-wide) mean rather than a group-specific mean
#' @param groupvar string representing the name of the variable denoting which sample belongs to which group.  Can be \code{"none"} (if you want the study-wide mean), or must correspond to the name of a column of \code{pData(gown)}.  Usually a categorical variable.
#' @param groupname string representing which group's expression means you want to plot.  Can be \code{"none"} (if you want the study-wide mean), \code{"all"} (if you want a multipanel plot of each group's mean expression), or any of the levels of \code{groupvar}.
#' @param colorby one of \code{"transcript"} or \code{"exon"}, indicating which feature's abundances should dictate plot coloring. 
#' @param legend if \code{TRUE} (as it is by default), a color legend is drawn on top of the plot indicating the scale for feature abundances.
#' @param labelTranscripts if \code{TRUE}, transcript ids are labeled on the left side of the plot. Default \code{FALSE}.
#' @return produces a plot of the transcript structure for the specified gene in the current graphics device, colored by study-wide or group-specific mean expression level.
#' @seealso \code{\link{plotTranscripts}} 
#' @author Alyssa Frazee
#' @export

plotMeans = function(gene, gown, 
    overall = TRUE,
    groupvar = "none", 
    groupname = "none", 
	dattype = c("cov", "FPKM", "rcount", "ucount", "mrcount", "mcov"),
	colorby = c("transcript", "exon"),
    legend = TRUE,
    labelTranscripts = FALSE){
    
    dattype = match.arg(dattype)
    colorby = match.arg(colorby)

    if(colorby == "transcript" & !(dattype %in% c("cov", "FPKM"))){
        stop("transcripts only have cov and FPKM measurements")
    }
    if(colorby == "exon" & dattype=="FPKM"){
        stop("exons do not have FPKM measurements")
    }
    if(class(gown)!="ballgown") stop("gown must be a ballgown object")

    if(!overall & (groupvar=="none"|groupname=="none")){
        stop("to plot means for a specific group, please provide both the grouping variable name and which group you're plotting for")
    }
    if(groupname == "all" & overall){
        warning("Plotting the study-wide mean for each feature. To plot each group's mean separately, set overall=FALSE.")
    }
    if(groupname == "all" & groupvar == "none"){
        stop("to plot means for all groups, please provide groupvar (the name of the group variable)")
    }
    
    # some setup:
    ma = as.data.frame(structure(gown)$trans)
    thetranscripts = indexes(gown)$t2g$t_id[indexes(gown)$t2g$g_id==gene]
    thetranscripts = paste0("tx", thetranscripts)
    gtrans = subset(ma, element %in% thetranscripts)
    gtrans$tid = as.numeric(sapply(gtrans$element, function(x) as.numeric(substr(x,3,nchar(x)))))
    xax = seq(min(gtrans$start), max(gtrans$end), by=1)
    numtx = length(unique(thetranscripts))
    if(groupvar != "none"){
        pdatacol = which(names(pData(gown))==groupvar) # the group variable (column index)
    }
    if(groupname == "all"){
        numplots = length(unique(pData(gown)[,pdatacol]))
        mf = c(ceiling(sqrt(numplots)), floor(sqrt(numplots)))
        plot_titles = unique(pData(gown)[,pdatacol])
    }else{
        numplots = 1
        mf = c(1,1)
        plot_titles = ifelse(groupname=="none", "mean across groups", groupname)
    }
    if(overall){
        samples = list(pData(gown)[,1])
    }else if(groupname == "all"){
        samples = split(pData(gown)[,1], pData(gown)[,pdatacol])
    }else{
        samples = list(pData(gown)[,1][pData(gown)[,pdatacol]==groupname]) 
    }
    
    # plot base:
    westval = ifelse(labelTranscripts, 4, 2)
    par(mar=c(5, westval, 4, 2), mfrow = mf)
    ymax = ifelse(legend, numtx+1.5, numtx+1)

    # the for-loop only has >1 iteration if groupname == "all"
    for(p in 1:numplots){
        plot(xax, rep(0,length(xax)), ylim=c(0,ymax), 
        type="n", xlab="genomic position", 
        main=paste0(gene,": ",plot_titles[p]), 
        yaxt = "n", ylab="")

        # calculate means
        if(colorby == "transcript"){
            # mean transcript-level expression measurement for this group
            smalldat = subset(texpr(gown), t_id %in% gtrans$tid)
            coltype = as.character(sapply(names(data(gown)$trans), gettype))
            colsamp = as.character(sapply(names(data(gown)$trans), getsamp))
            datacols = which(coltype==dattype & colsamp %in% samples[[p]])
            smalldat2 = smalldat[,datacols]
            tmeans = apply(smalldat2, 1, mean)
        }

        if(colorby == "exon"){
            # mean exon-level expression measurement for this group
            smalldat = subset(eexpr(gown), e_id %in% gtrans$id)
            coltype = as.character(sapply(names(data(gown)$exon), gettype))
            colsamp = as.character(sapply(names(data(gown)$exon), getsamp))
            datacols = which(coltype==dattype & colsamp %in% samples[[p]])
            smalldat2 = smalldat[,datacols]
            emeans = apply(smalldat2, 1, mean)
        }
    
        # make color scale:
        maxcol = quantile(as.matrix(smalldat[, which(coltype==dattype)]), 0.99)
        colscale = seq(0,maxcol,length.out=200)

        # draw the transcripts
        for(tx in unique(gtrans$tid)){
            if(colorby == "transcript") mycolor = closestColor(tmeans[which(smalldat$t_id==tx)], colscale)
            txind = which(unique(gtrans$tid)==tx)
            gtsub = gtrans[gtrans$tid==tx,]
            gtsub = gtsub[order(gtsub$start),]
            for(exind in 1:dim(gtsub)[1]){
                if(colorby == "exon") mycolor = closestColor(emeans[which(smalldat$e_id==gtsub$id[exind])], colscale)
                polygon(x=c(gtsub$start[exind], gtsub$start[exind], gtsub$end[exind], gtsub$end[exind]), y=c(txind-0.4,txind+0.4,txind+0.4,txind-0.4), col=mycolor)
                if(exind!=dim(gtsub)[1]) lines(c(gtsub$end[exind],gtsub$start[exind+1]),c(txind, txind), lty=2, col="gray60")
            }
        }

        # draw the legend:
        if(legend){
            leglocs = seq(min(xax)+1, max(xax)-1, length=length(colscale)+1)
            for(i in 1:length(colscale)){
                polygon(x = c(leglocs[i], leglocs[i], leglocs[i+1], leglocs[i+1]), y = c(ymax-0.3, ymax, ymax, ymax-0.3), border=NA, col = rev(heat.colors(length(colscale)))[i])
            }
            text(x = seq(min(xax)+1, max(xax)-1, length = 20), y = rep(ymax+0.1, 20), labels = round(colscale,2)[seq(1,length(colscale), length=20)], cex=0.5 )	
            text(x = median(xax), y = ymax-0.5, labels=paste("mean expression, by",  colorby), cex=0.5)
        }

        if(labelTranscripts) axis(side=2, at=c(1:numtx), labels=unique(gtrans$tid), cex.axis=0.75, las=1)
    }
}

