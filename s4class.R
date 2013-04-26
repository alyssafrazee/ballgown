## ballgown S4 class (other names welcome)

# class definition

setClass("ballgown", 
	representation(
	  data = "list",			# coverage data
	  indexes = "list",			# reference information
	  structure = "list",		# assembly information
	  dirs = "character",		# directories where data is stored
	  mergedDate = "character"	# date the object was created
	)
)

# constructor + helpers (this will be the slow-ish part)

## intron helper
.readIntron <- function(file){
  intron <- read.table(file, header=TRUE, sep="\t", colClasses=c("integer", "character", "factor", "integer", "integer", "integer", "integer", "numeric"))
  intron <- intron[order(intron$i_id), ]
  rownames(intron) <- 1:nrow(intron)
  return(intron)
}

## exon helper
.readExon <- function(file) {
  exon <- read.table(file, header=TRUE, sep="\t", colClasses=c("integer", "character", "factor", "integer", "integer", "integer", "integer", "numeric", "numeric", "numeric", "numeric", "numeric"))
  exon <- exon[order(exon$e_id), ]
  rownames(exon) <- 1:nrow(exon)
  return(exon)
}

## transcript helper
.readTrans <- function(file) {
  trans <- read.table(file, header=TRUE, sep="\t", colClasses=c("integer", "character", "factor", "integer", "integer", "character", "integer", "integer", "character", "character", "numeric", "numeric"))
  trans <- trans[order(trans$t_id), ]
  rownames(trans) <- 1:nrow(trans)
  return(trans)
}

ballgown = function(dataDir, samplePattern, bamfiles = NULL, pData = NULL, verbose=TRUE, ...) {
  if(verbose) message(date())

  ## Load required pkgs
  suppressMessages(library(plyr))
  suppressMessages(library(GenomicRanges))

  ## Identify the sample directories
  dirs <- list.files(path=dataDir, pattern=samplePattern, full.names=TRUE)
  names(dirs) <- list.files(path=dataDir, pattern=samplePattern)
  n <- length(dirs)

  ## Read linking tables
  if(verbose) message(paste0(date(), ": Reading linking tables"))
  e2t <- read.table(list.files(dirs[1], "e2t.ctab", full.names=TRUE), header=TRUE, sep="\t", colClasses=c("integer", "integer"))
  i2t <- read.table(list.files(dirs[1], "i2t.ctab", full.names=TRUE), header=TRUE, sep="\t", colClasses=c("integer", "integer"))

  ## Order by transcript id
  e2t <- e2t[order(e2t$t_id), ]
  i2t <- i2t[order(i2t$t_id), ]
  rownames(e2t) <- 1:nrow(e2t)
  rownames(i2t) <- 1:nrow(i2t)

  ## Read counts for all introns in <reference_transcripts>
  if(verbose) message(paste0(date(), ": Reading intron data files"))
  intronFiles <- sapply(dirs, list.files, pattern="i_data.ctab", full.names=TRUE)
  intronAll <- lapply(intronFiles, .readIntron)

  ## Merge the results
  if(verbose) message(paste0(date(), ": Merging intron data"))
  intron <- join_all(intronAll, by=c("i_id", "chr", "strand", "start", "end"), type="left")
  colnames(intron)  <- c("i_id", "chr", "strand", "start", "end", paste(c("rcount", "ucount", "mrcount"), rep(names(dirs), each=3), sep="."))

  ## Make intron data into GRanges object
  #A. fix strand information for compatibility w/ GRanges
  intron$strand = as.character(intron$strand)
  intron$strand[intron$strand=="."] <- "*"
  #B. get names of transcripts each intron belongs to
  tnamesin = split(i2t$t_id, i2t$i_id)
  #C. make the GRanges object
  introngr = GRanges(seqnames = Rle(intron$chr), ranges = IRanges(start=intron$start, end=intron$end), strand = Rle(intron$strand), id=intron$i_id, transcripts = as.character(tnamesin))

  ## Read read counts and raw coverage info for all exons in <reference_transcripts>
  if(verbose) message(paste0(date(), ": Reading exon data files"))
  exonFiles <- sapply(dirs, list.files, pattern="e_data.ctab", full.names=TRUE)
  exonAll <- lapply(exonFiles, .readExon)

  ## Read exon data
  if(verbose) message(paste0(date(), ": Merging exon data"))
  exon <- join_all(exonAll, by=c("e_id", "chr", "strand", "start", "end"), type="left")
  colnames(exon) <- c("e_id", "chr", "strand", "start", "end", paste(c("rcount", "ucount", "mrcount", "cov", "cov_sd", "mcov", "mcov_sd"), rep(names(dirs), each=7), sep="."))

  ## Make exon data into GRanges object
  #A. fix strand information for compatibility w/ GRanges
  exon$strand = as.character(exon$strand)
  exon$strand[exon$strand=="."] <- "*"
  #B. get names of transcripts each exon belongs to
  tnamesex = split(e2t$t_id, e2t$e_id)
  #C. make the GRanges object
  exongr = GRanges(seqnames = Rle(exon$chr), ranges = IRanges(start=exon$start, end=exon$end), strand = Rle(exon$strand), id=exon$e_id, transcripts = as.character(tnamesex))

  ## Read transcript data
  if(verbose) message(paste0(date(), ": Reading transcript data files"))
  transFiles <- sapply(dirs, list.files, pattern="t_data.ctab", full.names=TRUE)
  transAll <- lapply(transFiles, .readTrans)

  ## Merge the results
  if(verbose) message(paste0(date(),": Merging transcript data"))
  trans <- join_all(transAll, by=c("t_id", "chr", "strand", "start", "end", "t_name", "num_exons", "length", "gene_id", "gene_name"), type="left")
  colnames(trans) <- c("t_id", "chr", "strand", "start", "end", "t_name", "num_exons", "length", "gene_id", "gene_name", paste(c("cov", "FPKM"), rep(names(dirs), each=2), sep="."))

   ## Make transcripts into a GRanges list object
   transgrl = split(exongr[e2t$e_id], e2t$t_id)
   names(transgrl) = paste0("tx", names(transgrl))

   ## Connect transcripts to genes:
   t2g = data.frame(t_id = trans$t_id, g_id = trans$gene_id)

   ## Read phenotype table, if given:
   if(!is.null(pData)){
       if(verbose) message(paste0(date(),": Reading phenotype table"))
       phx = read.table(pData, stringsAsFactors=FALSE, ...)
       theorder = sapply(names(dirs), function(x) which(phx$dirname==x))
       phx = phx[theorder,]
   }
   if(is.null(pData)) phx = NULL

  if(verbose) message("Wrapping up the results")
  result <- new("ballgown", data = list(intron=intron, exon=exon, trans=trans), indexes=list(e2t=e2t, i2t=i2t, t2g=t2g, bamfiles = bamfiles, pData = phx), structure = list(intron = introngr, exon = exongr, trans = transgrl), dirs=dirs, mergedDate=date())

  if(verbose) message(date())
  return(result)
}


# defining the slot getters
setGeneric("structure", function(x) standardGeneric("structure"))
setMethod("structure", "ballgown", function(x) x@structure)
setGeneric("data", function(x) standardGeneric("data"))
setMethod("data", "ballgown", function(x) x@data)
setGeneric("indexes", function(x) standardGeneric("indexes"))
setMethod("indexes", "ballgown", function(x) x@indexes)
setGeneric("dirs", function(x) standardGeneric("dirs"))
setMethod("dirs", "ballgown", function(x) x@dirs)
setGeneric("mergedDate", function(x) standardGeneric("mergedDate"))
setMethod("mergedDate", "ballgown", function(x) x@mergedDate)


# the show method:
setMethod("show", "ballgown", 
	function(object)
		cat(class(object), "instance with", length(structure(object)$trans), "assembled transcripts")
)


# define slot setters (still need to finish)
setGeneric("indexes<-", function(x, value) standardGeneric("indexes<-"))
setReplaceMethod("indexes", "ballgown", function(x, value) {x@indexes <- value; x})
###### AFTER DEFINING A VALIDITY METHOD:
setReplaceMethod("indexes", "ballgown", function(x, value) {x@indexes <- value; validObject(x); x})


# define coercion methods (still need to do)


# here is a subset method - can subset by anything in the transcript table.
setMethod("subset", "ballgown", function(x, cond){
	ctext = deparse(substitute(cond))
	trans = subset(data(x)$trans, eval(parse(text=ctext)))
	
	thetx = trans$t_id
	
	inttmp = split(indexes(x)$i2t$i_id, indexes(x)$i2t$t_id)
	theint = as.numeric(unique(unlist(inttmp[names(inttmp) %in% thetx])))
	intron = subset(data(x)$intron, i_id %in% theint)
	
	extmp = split(indexes(x)$e2t$e_id, indexes(x)$e2t$t_id)
	theex = as.numeric(unique(unlist(extmp[names(extmp) %in% thetx])))
	exon = subset(data(x)$exon, e_id %in% theex)
	
	e2t = subset(indexes(x)$e2t, t_id %in% thetx)
	i2t = subset(indexes(x)$i2t, t_id %in% thetx)
	t2g = subset(indexes(x)$t2g, t_id %in% thetx)
	
	introngr = structure(x)$intron[elementMetadata(structure(x)$intron)$id %in% theint]
	exongr = structure(x)$exon[elementMetadata(structure(x)$exon)$id %in% theex]
	grltxids = sapply(names(structure(x)$trans), function(a) as.numeric(substr(a, 3, nchar(a))))
	transgrl = structure(x)$trans[grltxids %in% thetx]
	
	return(new("ballgown", data = list(intron=intron, exon=exon, trans=trans), indexes = list(e2t=e2t, i2t=i2t, t2g=t2g, bamfiles = indexes(x)$bamfiles, pData = indexes(x)$pData), structure=list(intron = introngr, exon=exongr, trans=transgrl), dirs = dirs(x), mergedDate=mergedDate(x)))
} )




# plotTranscripts method: (I think this actually just needs to be a function, since "plotTranscripts" is not a current method)

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
	
	if(class(gown)!="ballgown") stop("gown must be a ballgown object")
	
	suppressMessages(library(GenomicRanges))
	
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

		}
	}
	
	ma = IRanges::as.data.frame(structure(gown)$trans)
	thetranscripts = indexes(gown)$t2g$t_id[indexes(gown)$t2g$g_id==gene]
	thetranscripts = paste0("tx", thetranscripts)
	gtrans = subset(ma, element %in% thetranscripts)
	gtrans$tid = as.numeric(sapply(gtrans$element, function(x) as.numeric(substr(x,3,nchar(x)))))
	xax = seq(min(gtrans$start), max(gtrans$end), by=1)
    numtx = length(unique(thetranscripts))
    par(mar=c(5,2,4,2))
    ymax = ifelse(legend, numtx+1.5, numtx+1)
    
    if(length(unique(gtrans$seqnames)) > 1) stop("Your gene appears to span multiple chromosomes, which is interesting but also kind of annoying, R-wise.  Please choose another gene until additional functionality is added!")
    if(length(unique(gtrans$strand)) > 1) stop("Your gene appears to contain exons from both strands, which is potentially interesting but also kind of confusing, so please choose another gene until we figure this sucker out.")
    
    # plot base:
    plot(xax, rep(0,length(xax)), ylim=c(0,ymax), type="n", xlab="genomic position", main=paste0(gene,": ",sampname), yaxt = "n", ylab="", )
    
    
    # set color scale:
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
    
    
    # draw the transcripts
    introntypes = unique(as.character(sapply(names(data(gown)$intron)[-c(1:5)], gettype)))
    color.introns = ifelse(gettype(samp) %in% introntypes, TRUE, FALSE)
    for(tx in unique(gtrans$tid)){
    	if(colorby == "transcript"){
    		mycolor = closestColor(data(gown)$trans[,col][which(data(gown)$trans$t_id==tx)], colscale)
    	}
    	txind = which(unique(gtrans$tid)==tx)
    	gtsub = gtrans[gtrans$tid==tx,]
    	gtsub = gtsub[order(gtsub$start),]
    	for(exind in 1:dim(gtsub)[1]){
    		if(colorby == "exon"){
    			mycolor = closestColor(data(gown)$exon[,col][which(data(gown)$exon$e_id==gtsub$id[exind])], colscale)
    		}
			polygon(x=c(gtsub$start[exind], gtsub$start[exind], gtsub$end[exind], gtsub$end[exind]), y=c(txind-0.4,txind+0.4,txind+0.4,txind-0.4), col=mycolor)
			if(exind!=dim(gtsub)[1]){
				if(!color.introns){
					lines(c(gtsub$end[exind],gtsub$start[exind+1]),c(txind, txind), lty=2, col="gray60")
				}
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
}



## and here is the plotMeans function:
plotMeans = function(gene, gown, groupvar, groupname, dattype = "cov", legend = TRUE, colorby = "transcript"){
	
	if(class(gown)!="ballgown") stop("gown must be a ballgown object")
	
	suppressMessages(library(GenomicRanges))
	
	# some setup:
	outcomecol = which(names(indexes(gown)$pData)==groupvar)
	ma = IRanges::as.data.frame(structure(gown)$trans)
	thetranscripts = indexes(gown)$t2g$t_id[indexes(gown)$t2g$g_id==gene]
	thetranscripts = paste0("tx", thetranscripts)
	gtrans = subset(ma, element %in% thetranscripts)
	gtrans$tid = as.numeric(sapply(gtrans$element, function(x) as.numeric(substr(x,3,nchar(x)))))
	xax = seq(min(gtrans$start), max(gtrans$end), by=1)
    numtx = length(unique(thetranscripts))
    par(mar=c(5,2,4,2))
    ymax = ifelse(legend, numtx+1.5, numtx+1)
    pdatacol = which(names(indexes(gown)$pData)==groupvar) # the group variable (column index)
	samples = indexes(gown)$pData$dirname[indexes(gown)$pData[,pdatacol]==groupname] # names of the samples in that group

    
    # plot base:
    plot(xax, rep(0,length(xax)), ylim=c(0,ymax), type="n", xlab="genomic position", main=paste0(gene,": ",groupname), yaxt = "n", ylab="", )


	# calculate means
	if(colorby == "transcript"){
		# mean transcript-level expression measurement for this group
		smalldat = subset(data(gown)$trans, t_id %in% gtrans$tid)
		coltype = as.character(sapply(names(data(gown)$trans), gettype))
		colsamp = as.character(sapply(names(data(gown)$trans), getsamp))
		datacols = which(coltype==dattype & colsamp %in% samples)
		smalldat2 = smalldat[,datacols]
		tmeans = apply(smalldat2, 1, mean)
		}
	
    if(colorby == "exon"){
    	# mean exon-level expression measurement for this group
    	smalldat = subset(data(gown)$exon, e_id %in% gtrans$id)
    	coltype = as.character(sapply(names(data(gown)$exon), gettype))
    	colsamp = as.character(sapply(names(data(gown)$exon), getsamp))
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




setwd("~/Google Drive/hopkins/research/_cufflinks visualization project")
load("small_gown.rda")
tinygown = new("ballgown", data=gown$data, indexes=gown$indexes, structure=gown$structure, dirs=gown$dirs, mergedDate=gown$mergedDate)
# some fixeruppers (all changed in the new read function)
indexes(tinygown)$t2g$g_id[indexes(tinygown)$t2g$t_id %in% data(tinygown)$trans$t_id] <- data(tinygown)$trans$gene_id
indexes(tinygown)$pData <- read.table("~/Desktop/phenotypes_all.txt", stringsAsFactors=FALSE, sep="\t", header=TRUE)
indexes(tinygown)$pData$dirname[18] <- NA
theorder = sapply(names(dirs(tinygown)), function(x) which(indexes(tinygown)$pData$dirname==x))
indexes(tinygown)$pData = indexes(tinygown)$pData[theorder,]

plotTranscripts("XLOC_000002", "rcount.orbFrontalF2", tinygown, colorby="exon")
plotTranscripts("XLOC_000002", "FPKM.orbFrontalF2", tinygown, colorby="transcript")
plotTranscripts("XLOC_000011", "rcount.orbFrontalF2", tinygown, colorby="exon")
plotTranscripts("XLOC_000011", "cov.orbFrontalF2", tinygown, colorby="transcript")

plotMeans("XLOC_000002", tinygown, groupvar = "outcome", groupname = "bipolar", dattype = "cov", legend = TRUE, colorby = "exon")
plotMeans("XLOC_000002", tinygown, groupvar = "outcome", groupname = "control", dattype = "cov", legend = TRUE, colorby = "transcript")
plotMeans("XLOC_000002", tinygown, groupvar = "outcome", groupname = "schizophrenia", dattype = "cov", legend = TRUE, colorby = "transcript")
plotMeans("XLOC_000002", tinygown, groupvar = "outcome", groupname = "depression", dattype = "cov", legend = TRUE, colorby = "transcript")


######## MERGING TRANSCRIPTS #########
# (1) define function calculating transcript overlaps within a gene.


bpset = function(tra){
	unique(unlist(sapply(1:length(ranges(tra)), function(i) c(start(tra[i,]):end(tra[i,])))))
}

transcriptOverlaps = function(gene, gown){
	
	txnames = indexes(gown)$t2g$t_id[indexes(gown)$t2g$g_id == gene]
	strucnames = as.numeric(substr(names(structure(gown)$trans),3,nchar(names(structure(gown)$trans))))
	inds = which(strucnames %in% txnames)
	tx = structure(gown)$trans[inds]
	
	overlapMat = matrix(NA, nrow=length(inds), ncol=length(inds))
	rownames(overlapMat) = colnames(overlapMat) = names(tx)
	diag(overlapMat) = 1
	for(ii in 1:length(inds)){
		for(jj in 1:length(inds)){
			if(ii==jj) next
			bpsi = bpset(tx[[ii]])
			bpsj = bpset(tx[[jj]])
			overlapMat[ii,jj] = length(intersect(bpsi, bpsj))/length(bpsi)
		}
	}
	
	return(overlapMat)
}  #number in row i, column j answers the question "what percentage of transcript i is overlapped by transcript j?"


system.time(transcriptOverlaps("XLOC_000002", tinygown)) #8 seconds elapsed time
system.time(transcriptOverlaps("XLOC_000011", tinygown)) #15 seconds elapsed time
transcriptOverlaps("XLOC_000011", tinygown)
plotMeans("XLOC_000011", tinygown, "outcome", "control", colorby="exon")


#(2) decide which transcripts to combine into one.  (I vote for the "reduce" function but what do I know?)



########### DIFFERENTIAL EXPRESSION TESTS ###########
# here is a way we can do differential expression (the test formerly known as "cufffix")
difftest = function(gown, feature = "transcript", method = "DESeq", dattype = "cov", ...){
	if(feature=="transcript"){
		if(dattype!="cov" & dattype!="FPKM") stop("transcripts only have cov and FPKM measurements")
		coltypes = as.character(sapply(names(data(gown)$trans), gettype))
		inds = which(coltypes==dattype)
		tab = data(gown)$trans[,inds]
	}
	if(feature=="exon"){
		exontypes = unique(as.character(sapply(names(data(gown)$exon)[-c(1:5)], gettype)))
		if(!(dattype %in% exontypes)) stop(paste0("exons only have the following measurements: ", paste(exontypes, collapse=", ")))
		coltypes = as.character(sapply(names(data(gown)$exon), gettype))
		inds = which(coltypes==dattype)
		tab = data(gown)$exon[,inds]
	}
	if(feature=="junction" | feature=="intron"){
		introntypes = c("rcount", "ucount", "mrcount")
		if(!(dattype %in% introntypes)) stop(paste0("introns/junctions only have the following measurements: ", paste(introntypes, collapse=", ")))
		coltypes = as.character(sapply(names(data(gown)$intron), gettype))
		inds = which(coltypes==dattype)
		tab = data(gown)$intron[,inds]
	}
	if(!(feature %in% c("transcript","exon","junction","intron"))) stop(paste0("unknown feature: ",feature,". Please choose one of transcript, exon, intron, or junction.  (\"intron\" and \"junction\" are the same test)."))
	
	# now do stuff!
	
}
