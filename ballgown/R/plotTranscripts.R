plotTranscripts = function(gene, samp = NULL, gown, 
	legend = TRUE, labelTranscripts = FALSE, 
	colorby = c("transcript", "exon", "none"),
	main = NULL){

	if(class(gown)!="ballgown") stop("gown must be a ballgown object")
	if(colorby!="none" & is.null(samp)) stop("to color by transcript or exon abundances, you must provide a specific sample. (use names(data(gown)$trans) or names(data(gown)$exon) to see sample names).")
	
	colorby = match.arg(colorby)

	if(colorby=="none") legend = FALSE
	if(colorby!="none"){
	# if "samp" is the character name of the column:
		if(is.character(samp)){
			if(colorby == "transcript"){
				col = which(names(data(gown)$trans)==samp)
				if(gettype(samp)!="cov" & gettype(samp)!="FPKM") stop("transcripts only have cov and FPKM measurements")
			}
			if(colorby == "exon"){
				exontypes = unique(as.character(sapply(names(data(gown)$exon)[-c(1:5)], gettype)))
				if(!(gettype(samp) %in% exontypes)) stop(paste0("exons only have the following measurements: ", paste(exontypes, collapse=", ")))
				col = which(names(data(gown)$exon)==samp)
			}
			sampname = samp
		}
		# if it is the column index:
		if(is.numeric(samp)){
			col = samp
			if(colorby == "transcript"){
				sampname = names(data(gown)$trans)[samp]	
				if(gettype(sampname)!="cov" & gettype(sampname)!="FPKM") stop("transcripts only have cov and FPKM measurements")
			}
			if(colorby == "exon"){
				sampname = names(data(gown)$exon)[samp]
				if(!(gettype(sampname) %in% exontypes)) stop(paste0("exons only have the following measurements: ", paste(exontypes, collapse=", ")))
			} #end if colorby=="exon"
		} #end is.numeric(samp)
	} #end if colorby != "none"

	ma = IRanges::as.data.frame(structure(gown)$trans)
	thetranscripts = indexes(gown)$t2g$t_id[indexes(gown)$t2g$g_id==gene]
	thetranscripts = paste0("tx", thetranscripts)
	gtrans = subset(ma, element %in% thetranscripts)
	gtrans$tid = as.numeric(sapply(gtrans$element, function(x) as.numeric(substr(x,3,nchar(x)))))
	xax = seq(min(gtrans$start), max(gtrans$end), by=1)
	numtx = length(unique(thetranscripts))
	westval = ifelse(labelTranscripts, 4, 2)
	par(mar=c(5, westval, 4, 2))
	ymax = ifelse(legend, numtx+1.5, numtx+1)
	
	if(length(unique(gtrans$seqnames)) > 1) stop("Your gene appears to span multiple chromosomes, which is interesting but also kind of annoying, R-wise.  Please choose another gene until additional functionality is added!")
	if(length(unique(gtrans$strand)) > 1) stop("Your gene appears to contain exons from both strands, which is potentially interesting but also kind of confusing, so please choose another gene until we figure this sucker out.")
    
	# plot base:
	plot(xax, rep(0,length(xax)), ylim=c(0,ymax), type="n", xlab="genomic position", yaxt = "n", ylab="")
	if(!is.null(main)){
		mytitle = main
	}else{
		if(colorby!="none"){
			mytitle = paste0(gene,": ",sampname)
		}else{
			mytitle = gene
		}
	}
	title(mytitle)

	# set color scale
	if(colorby!="none"){
		sampcoltype = gettype(samp)
		if(colorby == "transcript"){
			smalldat = subset(data(gown)$trans, gene_id==gene)
			coltype = as.character(sapply(names(data(gown)$trans), gettype))
		}
		if(colorby == "exon"){
			smalldat = subset(data(gown)$exon, e_id %in% gtrans$id)
			coltype = as.character(sapply(names(data(gown)$exon), gettype))
		}
		maxcol = quantile(as.matrix(smalldat[,coltype==sampcoltype]), 0.99)
		colscale = seq(0,maxcol,length.out=200)
		introntypes = unique(as.character(sapply(names(data(gown)$intron)[-c(1:5)], gettype)))
		color.introns = ifelse(gettype(samp) %in% introntypes, TRUE, FALSE)
	}
	if(colorby=="none") color.introns = FALSE

	# draw transcripts
	for(tx in unique(gtrans$tid)){
		if(colorby == "transcript") mycolor = closestColor(data(gown)$trans[,col][which(data(gown)$trans$t_id==tx)], colscale)
		if(colorby == "none") mycolor = "gray70"
		txind = which(unique(gtrans$tid)==tx)
		gtsub = gtrans[gtrans$tid==tx,]
		gtsub = gtsub[order(gtsub$start),]
		for(exind in 1:dim(gtsub)[1]){
			if(colorby == "exon") mycolor = closestColor(data(gown)$exon[,col][which(data(gown)$exon$e_id==gtsub$id[exind])], colscale)
			polygon(x=c(gtsub$start[exind], gtsub$start[exind], gtsub$end[exind], gtsub$end[exind]), y=c(txind-0.4,txind+0.4,txind+0.4,txind-0.4), col=mycolor)
			if(exind!=dim(gtsub)[1]){
				if(!color.introns) lines(c(gtsub$end[exind],gtsub$start[exind+1]),c(txind, txind), lty=2, col="gray60")
				if(color.introns){
					intronindex = which(data(gown)$intron$start == gtsub$end[exind]+1 & data(gown)$intron$end == gtsub$start[exind+1]-1 & data(gown)$intron$chr==unique(gtsub$seqnames) & data(gown)$intron$strand == unique(gtsub$strand))
					icolumnind = which(names(data(gown)$intron) == samp)
					icol = closestColor(data(gown)$intron[intronindex,icolumnind], colscale)
					lines(c(gtsub$end[exind]+10,gtsub$start[exind+1]-10),c(txind, txind), lwd=3, col=icol)
					lines(c(gtsub$end[exind],gtsub$start[exind+1]),c(txind+0.1, txind+0.1), lwd=0.5, col="gray60")
					lines(c(gtsub$end[exind],gtsub$start[exind+1]),c(txind-0.1, txind-0.1), lwd=0.5, col="gray60")
				}
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

	# label the transcripts on the y-axis
	if(labelTranscripts){
		axis(side=2, at=c(1:numtx), labels=unique(gtrans$tid), cex.axis=0.75, las=1)
	}
}
