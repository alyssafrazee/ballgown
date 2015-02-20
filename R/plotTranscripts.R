#' visualize structure of assembled transcripts
#'
#' @param gene name of gene whose transcripts will be plotted.  When using 
#'   Cufflinks output, usually of the form \code{"XLOC_######"}
#' @param gown ballgown object containing experimental and phenotype data
#' @param samples vector of sample(s) to plot. Can be \code{'none'} if only one 
#'   plot (showing transcript structure in gray) is desired. Use 
#'   \code{sampleNames(gown)} to see sample names for \code{gown}. Defaults to 
#'   \code{sampleNames(gown)[1]}.
#' @param colorby one of \code{"transcript"}, \code{"exon"}, or \code{"none"}, 
#'   indicating which feature's abundances should dictate plot coloring.  If 
#'   \code{"none"}, all transcripts are drawn in gray.
#' @param meas which expression measurement to color features by, if any. Must 
#'   match an available measurement for whatever feature you're plotting.
#' @param legend if \code{TRUE} (as it is by default), a color legend is drawn
#'   on top of the plot indicating scales for feature abundances.
#' @param labelTranscripts if \code{TRUE}, transcript ids are labeled on the 
#'   left side of the plot. Default \code{FALSE}.
#' @param main optional string giving the desired plot title.
#' @param blackBorders if \code{TRUE}, exon borders are drawn in black. 
#'   Otherwise, they are drawn in the same color as their transcript or exon. 
#'   Switching blackBorders to FALSE can be useful for visualizing abundances 
#'   for skinny exons and/or smaller plots, which can be the case when 
#'   \code{length(samples)} is large.
#' @param log if \code{TRUE}, color transcripts on the log scale. Default 
#'   \code{FALSE}. To account for expression values of 0, we add 1 to all 
#'   expression values before taking the log.
#' @param logbase log base to use if \code{log = TRUE}. Default 2.
#' @param customCol an optional vector of custom colors to color transcripts by.
#'   There must be the same number of colors as transcripts in the gene being 
#'   plotted.
#' @param customOrder an optional vector of transcript ids (matching ids in 
#'   \code{texpr(gown, 'all')$t_id}), indicating which order transcripts will 
#'   appear in the plot. All transcripts in \code{gene} must appear in the 
#'   vector exactly once.
#' 
#' @return produces a plot of the transcript structure for the specified gene in
#'   the current graphics device.
#' @seealso \code{\link{plotMeans}}, \code{\link{plotLatentTranscripts}}
#' @author Alyssa Frazee
#' @export
#' @examples \donttest{
#' data(bg)
#' 
#' # plot one gene for one sample:
#' plotTranscripts(gene='XLOC_000454', gown=bg, samples='sample12', meas='FPKM',
#'     colorby='transcript', 
#'     main='transcripts from gene XLOC_000454: sample 12, FPKM')
#' 
#' # plot one gene for many samples:
#' plotTranscripts('XLOC_000454', bg, 
#'     samples=c('sample01', 'sample06', 'sample12', 'sample19'), 
#'     meas='FPKM', colorby='transcript')
#' 
#' }
plotTranscripts = function(gene, gown, samples=NULL, colorby='transcript',
    meas='FPKM', legend=TRUE, labelTranscripts=FALSE, main=NULL, 
    blackBorders=TRUE, log=FALSE, logbase=2, customCol=NULL, customOrder=NULL){

    if(class(gown)!="ballgown") stop("gown must be a ballgown object")

    if(gown@RSEM & colorby=='exon'){
        stop(.makepretty('RSEM objects do not yet have exon-level measurements,
            so must color by transcript.'))
    }

    if(!gown@RSEM & meas=='TPM'){
        stop('only RSEM objects have TPM measurements.')
    }

    if(is.null(samples)){
        samples = sampleNames(gown)[1]
        if(colorby!='none') message(paste('defaulting to sample', samples))
    }
    
    stopifnot(colorby %in% c('transcript', 'exon', 'none'))

    if(colorby == 'transcript'){
        stopifnot(meas %in% c('cov', 'FPKM', 'TPM'))
    }
    if(colorby == 'exon'){
        emeas = c('rcount', 'ucount', 'mrcount', 'cov', 'cov_sd', 
            'mcov', 'mcov_sd')
        stopifnot(meas %in% emeas)
    }

    if(!(gown@meas == 'all' | meas %in% gown@meas)){
        stop(paste('gown does not contain', meas, 'measurements.'))
    }

    if(colorby=="none") legend = FALSE
    
    if(!is.null(customCol) & (colorby!="transcript")){
        stop("Custom coloring is only available at transcript level currently")
    }

    if(!is.null(customCol) & legend){
        stop('legend must be FALSE if you provide custom colors')
    }

    n = length(samples)
    westval = ifelse(labelTranscripts, 4, 2)
    if(n > 1){
        numrows = floor(sqrt(n))
        numcols = ceiling(n/numrows)
        par(mfrow=c(numrows,numcols), mar=c(5,westval,4,2), oma=c(0, 0, 2, 0))
    }else{
        par(mar=c(5, westval, 4, 2))
    }

    ma = IRanges::as.data.frame(structure(gown)$trans)
    if(names(ma)[2] != 'group_name'){
        stop('IRanges::as.data.frame has changed. Please report as issue.')
    }
    thetranscripts = indexes(gown)$t2g$t_id[indexes(gown)$t2g$g_id==gene]
    if(!is.null(customOrder)){
        if(!all(sort(customOrder) == sort(thetranscripts))){
            stop(.makepretty('customOrder must include each transcript in gene
                exactly once.'))
        }
    }
    
    if(substr(ma$group_name[1],1,2) == "tx"){
        warning('your ballgown object was built with a deprecated version of
            ballgown - would probably be good to re-build!')
        thetranscripts = paste0('tx', thetranscripts)
    }
    
    if(!is.null(customCol) & (length(customCol)!=length(thetranscripts))){
        stop("You must have the same number of custom colors as transcripts")
    }

    gtrans = ma[ma$group_name %in% thetranscripts,]
    xax = seq(min(gtrans$start), max(gtrans$end), by=1)
    numtx = length(unique(thetranscripts))
    ymax = ifelse(legend, numtx+1.5, numtx+1)
    
    if(length(unique(gtrans$seqnames)) > 1){ 
        stop("gene spans multiple chromosomes (?)")
    }
    if(length(unique(gtrans$strand)) > 1){ 
        stop("gene contains exons from both strands (?)")
    }

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
            colnames(smalldat) = snames
        }
        if(log){
            smalldat = log(smalldat+1, base=logbase)
        }
        maxcol = quantile(as.matrix(smalldat), 0.99)
        colscale = seq(0, maxcol, length.out=200)
        introntypes = c('ucount', 'rcount', 'mrcount')
        color.introns = meas %in% introntypes & 
            gown@meas %in% c('all', introntypes)
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
        if(colorby != 'none'){
            colIndex = which(colnames(smalldat) == colName)
        }

        # draw transcripts
        if(is.null(customOrder)){
            transcript_loop = unique(gtrans$group_name)
        }else{
            transcript_loop = customOrder
        }
        for(tx in transcript_loop){
            txind = which(transcript_loop==tx)
            if(colorby == 'transcript'){
                if(!is.null(customCol)){
                    mycolor = customCol[txind]
                }else{
                    mycolor = closestColor(smalldat[,colIndex][which(t_id==tx)],
                        colscale)    
                }
                stopifnot(length(mycolor) > 0)
            }else if(colorby == 'none'){
                mycolor = 'gray70'
            }
            
            gtsub = gtrans[gtrans$group_name==tx,]
            gtsub = gtsub[order(gtsub$start),]
            for(exind in 1:dim(gtsub)[1]){
                if(colorby == "exon"){ 
                    mycolor = closestColor(
                        smalldat[,colIndex][which(e_id==gtsub$id[exind])], 
                        colscale)
                    stopifnot(length(mycolor) > 0)
                }
                borderCol = ifelse(blackBorders, 'black', mycolor)
                polygon(x=c(gtsub$start[exind], gtsub$start[exind], 
                        gtsub$end[exind], gtsub$end[exind]), 
                        y=c(txind-0.4,txind+0.4,txind+0.4,txind-0.4), 
                        col=mycolor, border=borderCol)
                if(exind!=dim(gtsub)[1]){
                    if(!color.introns){ 
                        lines(c(gtsub$end[exind],gtsub$start[exind+1]), 
                            c(txind, txind), lty=2, col="gray60")
                    }
                    if(color.introns){
                        intronindex = which(
                            expr(gown)$intron$start == gtsub$end[exind]+1 & 
                            expr(gown)$intron$end == gtsub$start[exind+1]-1 & 
                            expr(gown)$intron$chr==unique(gtsub$seqnames) & 
                            expr(gown)$intron$strand == unique(gtsub$strand))
                        icolumnind = which(names(expr(gown)$intron) == colName)
                        icol = closestColor(
                            expr(gown)$intron[intronindex,icolumnind], colscale)
                        lines(c(gtsub$end[exind]+10,gtsub$start[exind+1]-10), 
                            c(txind, txind), lwd=3, col=icol)
                        lines(c(gtsub$end[exind],gtsub$start[exind+1]),
                            c(txind+0.1, txind+0.1), lwd=0.5, col="gray60")
                        lines(c(gtsub$end[exind],gtsub$start[exind+1]),
                            c(txind-0.1, txind-0.1), lwd=0.5, col="gray60")
                    }#end if color.introns
                }# end if exind != last exon
            }# end loop over exons
        }# end loop over transcripts
    
        # draw the legend:
        if(legend){
            leglocs = seq(min(xax)+1, max(xax)-1, length=length(colscale)+1)
            for(i in 1:length(colscale)){
                polygon(x=c(leglocs[i], leglocs[i], leglocs[i+1], leglocs[i+1]),
                    y=c(ymax-0.3, ymax, ymax, ymax-0.3), border=NA, 
                    col=rev(heat.colors(length(colscale)))[i])
            }
            text(x = seq(min(xax)+1, max(xax)-1, length=10), 
                y=rep(ymax+0.1, 10), 
                labels=round(colscale,2)[seq(1,length(colscale), length=10)], 
                cex=0.5) 
            if(log){
                text(x=median(xax), y=ymax-0.5, 
                    labels=paste("log expression, by",  colorby), cex=0.5)
            }else{
                text(x=median(xax), y=ymax-0.5, 
                    labels=paste("expression, by",  colorby), cex=0.5)
            }
        }

        # label the transcripts on the y-axis (if asked for)
        if(labelTranscripts){
            axis(side=2, at=c(1:numtx), labels=transcript_loop, cex.axis=0.75, 
                las=1)
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

