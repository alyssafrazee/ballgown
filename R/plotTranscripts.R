#' visualize structure of assembled transcripts
#'
#' @param gene name of gene whose transcripts will be plotted.  When using Cufflinks output, usually of the form \code{"XLOC_######"}
#' @param gown ballgown object containing experimental and phenotype data
#' @param samples vector of sample(s) to plot. Can be \code{'none'} if only one plot (showing transcript structure in gray) is desired. Use \code{sampleNames(gown)} to see sample names for \code{gown}. Defaults to \code{sampleNames(gown)[1]}.
#' @param colorby one of \code{"transcript"}, \code{"exon"}, or \code{"none"}, indicating which feature's abundances should dictate plot coloring.  If \code{"none"}, all transcripts are drawn in gray.
#' @param meas which expression measurement to color features by, if any. Must match an available measurement for whatever feature you're plotting.
#' @param legend if \code{TRUE} (as it is by default), a color legend is drawn on top of the plot indicating scales for feature abundances.
#' @param labelTranscripts if \code{TRUE}, transcript ids are labeled on the left side of the plot. Default \code{FALSE}.
#' @param main optional string giving the desired plot title.
#' @param colorBorders if \code{TRUE}, exon borders are also drawn in color (instead of black, as they are by default). Useful for visualizing abundances for skinny exons and/or smaller plots, as often happens when \code{length(samples)} is large.
#' @param log if \code{TRUE}, color transcripts on the log scale. Default \code{FALSE}. To account for expression values of 0, we add 1 to all expression values before taking the log.
#' @param logbase log base to use if \code{log = TRUE}. default 2.
#' @return produces a plot of the transcript structure for the specified gene in the current graphics device.
#' @seealso \code{\link{plotMeans}} 
#' @author Alyssa Frazee
#' @export

plotTranscripts = function(gene, gown, samples = NULL, 
    colorby = 'transcript', meas = 'FPKM', legend = TRUE, 
    labelTranscripts = FALSE, main = NULL, colorBorders = FALSE,
    log = FALSE, logbase = 2){


    if(class(gown)!="ballgown") stop("gown must be a ballgown object")
    if(is.null(samples)){
        samples = sampleNames(gown)[1]
        if(colorby!='none') message(paste('defaulting to sample',samples))
    }
    
    stopifnot(colorby %in% c('transcript', 'exon', 'none'))
    if(colorby == 'transcript'){
        stopifnot(meas %in% c('cov', 'FPKM'))
    }
    if(colorby == 'exon'){
        stopifnot(meas %in% c('rcount', 'ucount', 'mrcount', 'cov', 'cov_sd', 'mcov', 'mcov_sd'))
    }

    if(colorby=="none") legend = FALSE

    n = length(samples)
    westval = ifelse(labelTranscripts, 4, 2)
    if(n > 1){
        numrows = round(sqrt(n))
        numcols = ceiling(n/numrows)
        par(mfrow=c(numrows, numcols), mar=c(5, westval, 4, 2), oma = c(0, 0, 2, 0))
    }else{
        par(mar=c(5, westval, 4, 2))
    }

    ma = IRanges::as.data.frame(structure(gown)$trans)
    thetranscripts = indexes(gown)$t2g$t_id[indexes(gown)$t2g$g_id==gene]
    if(substr(ma$element[1] == "tx"){
        warning('your ballgown object was built with a deprecated version of ballgown - would probably be good to re-build!')
        thetranscripts = paste0('tx',thetranscripts)
    }
    gtrans = subset(ma, element %in% thetranscripts)
    xax = seq(min(gtrans$start), max(gtrans$end), by=1)
    numtx = length(unique(thetranscripts))
    ymax = ifelse(legend, numtx+1.5, numtx+1)
    
    if(length(unique(gtrans$seqnames)) > 1) stop("Your gene appears to span multiple chromosomes, which is interesting but also kind of annoying, R-wise.  Please choose another gene until additional functionality is added!")
    if(length(unique(gtrans$strand)) > 1) stop("Your gene appears to contain exons from both strands, which is potentially interesting but also kind of confusing, so please choose another gene until we figure this sucker out.")

    # set appropriate color scale:
    if(colorby != 'none'){
        g_id = texpr(gown, 'all')$gene_id
        if(colorby == "transcript"){
            smalldat = texpr(gown, meas)[which(g_id == gene),]
            t_id = texpr(gown, 'all')$t_id[which(g_id == gene)]
        }
        if(colorby == "exon"){
            e_id_full = eexpr(gown, 'all')$e_id
            smalldat = eexpr(gown, meas)[which(e_id_full %in% gtrans$id),]
            e_id = e_id_full[which(e_id_full %in% gtrans$id)]
        }
        if(numtx == 1){
            snames = names(smalldat)
            smalldat = matrix(smalldat, nrow=1)
            names(smalldat) = snames
        }
        if(log){
            smalldat = log(smalldat+1, base=logbase)
        }
        maxcol = quantile(as.matrix(smalldat), 0.99)
        colscale = seq(0, maxcol, length.out=200)
        introntypes = unique(as.character(sapply(names(data(gown)$intron)[-c(1:5)], gettype)))
        color.introns = meas %in% introntypes
    }else{
        color.introns = FALSE
    }

    # make plot(s) (one for each sample):
    for(s in 1:n){
        ## make base plot:
        plot(xax, rep(0,length(xax)), ylim=c(0,ymax), 
            type="n", xlab="genomic position", yaxt = "n", ylab="")
        if(n > 1){
            title(samples[s])
        }
        
        colName = paste(meas, samples[s], sep='.')

        # draw transcripts
        for(tx in unique(gtrans$element)){
            if(colorby == "transcript"){
                colIndex = which(names(smalldat) == colName)
                mycolor = closestColor(smalldat[,colIndex][which(t_id==tx)], colscale)
            }
            if(colorby == "none") mycolor = "gray70"
            txind = which(unique(gtrans$element)==tx)
            gtsub = gtrans[gtrans$element==tx,]
            gtsub = gtsub[order(gtsub$start),]
            for(exind in 1:dim(gtsub)[1]){
                if(colorby == "exon") mycolor = closestColor(smalldat[,colIndex][which(e_id==gtsub$id[exind])], colscale)
                borderCol = ifelse(colorBorders, mycolor, 'black')
                polygon(x=c(gtsub$start[exind], gtsub$start[exind], 
                    gtsub$end[exind], gtsub$end[exind]), 
                    y=c(txind-0.4,txind+0.4,txind+0.4,txind-0.4), 
                    col=mycolor, border=borderCol)
                if(exind!=dim(gtsub)[1]){
                    if(!color.introns) lines(c(gtsub$end[exind],gtsub$start[exind+1]),c(txind, txind), lty=2, col="gray60")
                    if(color.introns){
                        intronindex = which(data(gown)$intron$start == gtsub$end[exind]+1 & data(gown)$intron$end == gtsub$start[exind+1]-1 & data(gown)$intron$chr==unique(gtsub$seqnames) & data(gown)$intron$strand == unique(gtsub$strand))
                        icolumnind = which(names(data(gown)$intron) == colName)
                        icol = closestColor(data(gown)$intron[intronindex,icolumnind], colscale)
                        lines(c(gtsub$end[exind]+10,gtsub$start[exind+1]-10),c(txind, txind), lwd=3, col=icol)
                        lines(c(gtsub$end[exind],gtsub$start[exind+1]),c(txind+0.1, txind+0.1), lwd=0.5, col="gray60")
                        lines(c(gtsub$end[exind],gtsub$start[exind+1]),c(txind-0.1, txind-0.1), lwd=0.5, col="gray60")
                    }#end if color.introns
                }# end if exind != last exon
            }# end loop over exons
        }# end loop over transcripts
    
        # draw the legend:
        if(legend){
            leglocs = seq(min(xax)+1, max(xax)-1, length=length(colscale)+1)
            for(i in 1:length(colscale)){
                polygon(x = c(leglocs[i], leglocs[i], leglocs[i+1], leglocs[i+1]), y = c(ymax-0.3, ymax, ymax, ymax-0.3), border=NA, col = rev(heat.colors(length(colscale)))[i])
            }
            text(x = seq(min(xax)+1, max(xax)-1, length = 10), y = rep(ymax+0.1, 10), labels = round(colscale,2)[seq(1,length(colscale), length=10)], cex=0.5 ) 
            if(log){
                text(x = median(xax), y = ymax-0.5, labels=paste("log expression, by",  colorby), cex=0.5)
            }else{
                text(x = median(xax), y = ymax-0.5, labels=paste("expression, by",  colorby), cex=0.5)                
            }
        }

        # label the transcripts on the y-axis (if asked for)
        if(labelTranscripts){
            axis(side=2, at=c(1:numtx), labels=unique(gtrans$element), cex.axis=0.75, las=1)
        }

    } #end loop over samples

    # put title on plot
    if(!is.null(main)){
        title(main, outer=(n>1))
    }else{
        if(n==1){
            title(paste0(gene,': ',samples))
        }else{
            title(gene, outer=TRUE)
        }
    }

}

