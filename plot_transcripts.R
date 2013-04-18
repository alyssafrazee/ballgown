## visualizing ballgown/cufflinks output
## AF 11 April 2013

setwd("~/Google Drive/hopkins/research/_cufflinks visualization project")
## load a smaller test dataset (representing full ballgown output)
load("small_gown.rda")  # made with "readGown.R" on cluster

# clean up (this has been fixed in cluster script)
gown$indexes$t2g$g_id[gown$indexes$t2g$t_id %in% gown$data$trans$t_id] = gown$data$trans$gene_id



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
plotTranscripts = function(gene, samp, gown, legend = TRUE, colorby = "transcript"){
	
	suppressMessages(library(GenomicRanges))
	
	# if "samp" is the character name of the column:
	if(is.character(samp)){
		if(colorby == "transcript"){
			col = which(names(gown$data$trans)==samp)
		}
		if(colorby == "exon"){
			col = which(names(gown$data$exon)==samp)
		}
		sampname = samp
	}
	# if it is the column index:
	if(is.numeric(samp)){
		col = samp
		if(colorby == "transcript"){
			sampname = names(gown$data$trans)[samp]	
		}
		if(colorby == "exon"){
			sampname = names(gown$data$exon)[samp]
		}
	}
	
	ma = IRanges::as.data.frame(gown$structure$trans)
	thetranscripts = gown$indexes$t2g$t_id[gown$indexes$t2g$g_id==gene]
	thetranscripts = paste0("tx", thetranscripts)
	gtrans = subset(ma, element %in% thetranscripts)
	gtrans$tid = as.numeric(sapply(gtrans$element, function(x) as.numeric(substr(x,3,nchar(x)))))
	xax = seq(min(gtrans$start), max(gtrans$end), by=1)
    numtx = length(unique(thetranscripts))
    par(mar=c(5,2,4,2))
    ymax = ifelse(legend, numtx+1.5, numtx+1)
    
    
    # plot base:
    plot(xax, rep(0,length(xax)), ylim=c(0,ymax), type="n", xlab="genomic position", main=paste0(gene,": ",sampname), yaxt = "n", ylab="", )
    
    
    # set color scale:
	sampcoltype = gettype(samp)
    if(colorby == "transcript"){
    	smalldat = subset(gown$data$trans, gene_id==gene)
		coltype = as.character(sapply(names(gown$data$trans), gettype))
    }
    if(colorby == "exon"){
	   	smalldat = subset(gown$data$exon, e_id %in% gtrans$id)
    	coltype = as.character(sapply(names(gown$data$exon), gettype))
    }
    maxcol = quantile(as.matrix(smalldat[,coltype==sampcoltype]), 0.99)
    colscale = seq(0,maxcol,length.out=200)
    
    
    # draw the transcripts
    for(tx in unique(gtrans$tid)){
    	if(colorby == "transcript"){
    		mycolor = closestColor(gown$data$trans[,col][which(gown$data$trans$t_id==tx)], colscale)
    	}
    	txind = which(unique(gtrans$tid)==tx)
    	gtsub = gtrans[gtrans$tid==tx,]
    	gtsub = gtsub[order(gtsub$start),]
    	for(exind in 1:dim(gtsub)[1]){
    		if(colorby == "exon"){
    			mycolor = closestColor(gown$data$exon[,col][which(gown$data$exon$e_id==gtsub$id[exind])], colscale)
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
		text(x = median(xax), y = ymax-0.5, labels=paste("expression, by",  colorby), cex=0.5)
		}
}



plotTranscripts("XLOC_000002", "cov.orbFrontalF2", gown, colorby="exon")
plotTranscripts("XLOC_000002", "cov.orbFrontalF2", gown, colorby="transcript")
plotTranscripts("XLOC_000011", "cov.orbFrontalF2", gown, colorby="exon")
plotTranscripts("XLOC_000011", "cov.orbFrontalF2", gown, colorby="transcript")

# pData wasn't loaded before, but re-running readGown.R will do it. I've loaded it locally:
gown$indexes$pData <- read.table("~/Desktop/phenotypes_all.txt", stringsAsFactors=FALSE, sep="\t", header=TRUE)
gown$indexes$pData$dirname[18] <- NA
theorder = sapply(names(gown$dirs), function(x) which(gown$indexes$pData$dirname==x))
gown$indexes$pData = gown$indexes$pData[theorder,]









