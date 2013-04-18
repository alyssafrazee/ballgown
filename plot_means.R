### plotMeans - allows user to color transcripts or exon by group means


### helper functions 
last = function(x) return(tail(x,n=1))

closestColor = function(x, colscale){
	choices = rev(heat.colors(length(colscale)))
	diffs = abs(x-colscale)
	return(choices[which.min(diffs)])
}

gettype = function(x) strsplit(x, split="\\.")[[1]][1]
getsamp = function(x) last(strsplit(x, split="\\.")[[1]])


### MAIN FUNCTION
plotMeans = function(gene, gown, groupvar, groupname, dattype = "cov", legend = TRUE, colorby = "transcript"){
	
	suppressMessages(library(GenomicRanges))
	
	# some setup:
	outcomecol = which(names(gown$indexes$pData)==groupvar)
	ma = IRanges::as.data.frame(gown$structure$trans)
	thetranscripts = gown$indexes$t2g$t_id[gown$indexes$t2g$g_id==gene]
	thetranscripts = paste0("tx", thetranscripts)
	gtrans = subset(ma, element %in% thetranscripts)
	gtrans$tid = as.numeric(sapply(gtrans$element, function(x) as.numeric(substr(x,3,nchar(x)))))
	xax = seq(min(gtrans$start), max(gtrans$end), by=1)
    numtx = length(unique(thetranscripts))
    par(mar=c(5,2,4,2))
    ymax = ifelse(legend, numtx+1.5, numtx+1)
    pdatacol = which(names(gown$indexes$pData)==groupvar) # the group variable (column index)
	samples = gown$indexes$pData$dirname[gown$indexes$pData[,pdatacol]==groupname] # names of the samples in that group

    
    # plot base:
    plot(xax, rep(0,length(xax)), ylim=c(0,ymax), type="n", xlab="genomic position", main=paste0(gene,": ",numtx," transcripts"), yaxt = "n", ylab="", )


	# calculate means
	if(colorby == "transcript"){
		# mean transcript-level expression measurement for this group
		smalldat = subset(gown$data$trans, t_id %in% gtrans$tid)
		coltype = as.character(sapply(names(gown$data$trans), gettype))
		colsamp = as.character(sapply(names(gown$data$trans), getsamp))
		datacols = which(coltype==dattype & colsamp %in% samples)
		smalldat2 = smalldat[,datacols]
		tmeans = apply(smalldat2, 1, mean)
		}
	
    if(colorby == "exon"){
    	# mean exon-level expression measurement for this group
    	smalldat = subset(gown$data$exon, e_id %in% gtrans$id)
    	coltype = as.character(sapply(names(gown$data$exon), gettype))
    	colsamp = as.character(sapply(names(gown$data$exon), getsamp))
		datacols = which(coltype==dattype & colsamp %in% samples)
		smalldat2 = smalldat[,datacols]
		emeans = apply(smalldat2, 1, mean)
        }
    
    # make color scale:
    maxcol = quantile(as.matrix(smalldat[, which(coltype==dattype)]), 0.99)
    colscale = seq(0,maxcol,length.out=200)


    # draw the transcripts
    for(tx in unique(gtrans$tid)){
    	if(colorby == "transcript"){
    		mycolor = closestColor(tmeans[which(smalldat$t_id==tx)], colscale)
    	}
    	txind = which(unique(gtrans$tid)==tx)
    	gtsub = gtrans[gtrans$tid==tx,]
    	gtsub = gtsub[order(gtsub$start),]
    	for(exind in 1:dim(gtsub)[1]){
    		if(colorby == "exon"){
    			mycolor = closestColor(emeans[which(smalldat$e_id==gtsub$id[exind])], colscale)
    		}
			polygon(x=c(gtsub$start[exind], gtsub$start[exind], gtsub$end[exind], gtsub$end[exind]), y=c(txind-0.4,txind+0.4,txind+0.4,txind-0.4), col=mycolor)
			if(exind!=dim(gtsub)[1]){
				lines(c(gtsub$end[exind],gtsub$start[exind+1]),c(txind, txind), lty=2, col="gray60")
			}	
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

}


plotMeans("XLOC_000002", gown, groupvar = "outcome", groupname = "bipolar", dattype = "cov", legend = TRUE, colorby = "exon")
plotMeans("XLOC_000002", gown, groupvar = "outcome", groupname = "control", dattype = "cov", legend = TRUE, colorby = "transcript")
plotMeans("XLOC_000002", gown, groupvar = "outcome", groupname = "schizophrenia", dattype = "cov", legend = TRUE, colorby = "transcript")
plotMeans("XLOC_000002", gown, groupvar = "outcome", groupname = "depression", dattype = "cov", legend = TRUE, colorby = "transcript")


