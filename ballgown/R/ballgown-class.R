### the ballgown class

setClass("ballgown", 
	representation(
	  data = "list",			# coverage data
	  indexes = "list",			# reference information
	  structure = "list",		# assembly information
	  dirs = "character",		# directories where ballgown data is stored
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

ballgown = function(dirs=NULL, dataDir=NULL, samplePattern=NULL, bamfiles = NULL, pData = NULL, verbose=TRUE, ...) {
  if(verbose) message(date())

  if(all(c(is.null(dirs),is.null(dataDir),is.null(samplePattern)))) stop("must provide either dirs or both dataDir and samplePattern")

  ## Determine where data is located
  if(is.null(dirs)){
    if(is.null(samplePattern)|is.null(dataDir)) stop("must provide one of dataDir or samplePattern if dirs is NULL.")
    dirs <- list.files(path=dataDir, pattern=samplePattern, full.names=TRUE)
    names(dirs) <- list.files(path=dataDir, pattern=samplePattern)
  }else{
    names(dirs) = sapply(dirs, function(x){
      tail(strsplit(x, split="/")[[1]],n=1)
    }, USE.NAMES=FALSE)
  }

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
  tnamesin.ord = as.character(tnamesin)[match(intron$i_id, names(tnamesin))]
  #C. make the GRanges object
  introngr = GRanges(seqnames = Rle(intron$chr), ranges = IRanges(start=intron$start, end=intron$end), strand = Rle(intron$strand), id=intron$i_id, transcripts = tnamesin.ord)

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
  tnamesex.ord = as.character(tnamesex)[match(exon$e_id, names(tnamesex))]
  exongr = GRanges(seqnames = Rle(exon$chr), ranges = IRanges(start=exon$start, end=exon$end), strand = Rle(exon$strand), id=exon$e_id, transcripts = tnamesex.ord)

  ## Read transcript data
  if(verbose) message(paste0(date(), ": Reading transcript data files"))
  transFiles <- sapply(dirs, list.files, pattern="t_data.ctab", full.names=TRUE)
  transAll <- lapply(transFiles, .readTrans)

  ## Merge the results
  if(verbose) message(paste0(date(),": Merging transcript data"))
  trans <- join_all(transAll, by=c("t_id", "chr", "strand", "start", "end", "t_name", "num_exons", "length", "gene_id", "gene_name"), type="left")
  colnames(trans) <- c("t_id", "chr", "strand", "start", "end", "t_name", "num_exons", "length", "gene_id", "gene_name", paste(c("cov", "FPKM"), rep(names(dirs), each=2), sep="."))

  ## Make transcripts into a GRanges list object
  mm = match(e2t$e_id, mcols(exongr)$id)
  transgrl = split(exongr[mm], e2t$t_id)
  names(transgrl) = paste0("tx", names(transgrl))

  ## Connect transcripts to genes:
  t2g = data.frame(t_id = trans$t_id, g_id = trans$gene_id)

  ## Read phenotype table, if given:
  if(is.character(pData)){
      if(verbose) message(paste0(date(),": Reading phenotype table"))
      phx = read.table(pData, stringsAsFactors=FALSE, ...)
      theorder = sapply(names(dirs), function(x) which(phx$dirname==x))
      phx = phx[theorder,]
  }
  if(is.data.frame(pData)){
       phx = pData
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
setGeneric("data<-", function(x, value) standardGeneric("data<-"))
setReplaceMethod("data", "ballgown", function(x, value) {x@data <- value; x})


###### AFTER DEFINING A VALIDITY METHOD:
#setReplaceMethod("indexes", "ballgown", function(x, value) {x@indexes <- value; validObject(x); x})


# define coercion methods (still need to do)


# here is a subset method - can subset by anything in the transcript table.
setGeneric("subset", function(x) standardGeneric("subset"))
setMethod("subset", "ballgown", function(x, cond, global=TRUE){

	ctext = ifelse(global, deparse(substitute(cond)), cond) # this means that inside another function, you can make a string argument to give to subset.
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


# methods for nice ways of accessing things (inspired by ExpressionSet syntax)
setGeneric("pData", function(x) standardGeneric("pData"))
setMethod("pData", "ballgown", function(x){
  return(indexes(x)$pData)
})
setGeneric("pData<-", function(x, value) standardGeneric("pData<-"))
setMethod("pData<-", "ballgown", function(x, value) {x@indexes$pData <- value; x})


setGeneric("texpr", function(x, ...) standardGeneric("texpr"))
setMethod("texpr", "ballgown", function(x, meas="all"){
  meas = match.arg(meas, c("cov","FPKM","all"))
  if(meas!="all"){
    expr = data(x)$trans[,-c(1:10)]
    expr = expr[,sapply(colnames(expr), function(x) strsplit(x,split="\\.")[[1]][1]==meas)]
    rownames(expr) = data(x)$trans$t_id
  }else{
    expr = data(x)$trans
  }
  return(expr)
})

setGeneric("eexpr", function(x, ...) standardGeneric("eexpr"))
setMethod("eexpr", "ballgown", function(x, meas="all"){
  meas = match.arg(meas, c("rcount","ucount","mrcount","cov","mcov","all"))
  if(meas!="all"){
    expr = data(x)$exon[,-c(1:5)]
    expr = expr[,sapply(colnames(expr), function(x) strsplit(x,split="\\.")[[1]][1]==meas)]
    rownames(expr) = data(x)$exon$e_id
  }else{
    expr = data(x)$exon
  }
  return(expr)
})

setGeneric("iexpr", function(x, ...) standardGeneric("iexpr"))
setMethod("iexpr", "ballgown", function(x, meas="all"){
  meas = match.arg(meas, c("rcount","ucount","mrcount","all"))
  if(meas!="all"){
    expr = data(x)$intron[,-c(1:5)]
    expr = expr[,sapply(colnames(expr), function(x) strsplit(x,split="\\.")[[1]][1]==meas)]
    rownames(expr) = data(x)$intron$i_id
  }else{
    expr = data(x)$intron
  }
  return(expr)
})

setGeneric("gexpr", function(x) standardGeneric("gexpr"))
setMethod("gexpr", "ballgown", function(x){
    gnames = indexes(x)$t2g$g_id
    inds_by_gene = split(seq(along=gnames), gnames)
    tmeas = texpr(x, "FPKM")
    gid_by_exon = lapply(1:nrow(texpr(x)), function(i){rep(texpr(x)$gene_id[i], texpr(x)$num_exons[i])})
    ulstruct = unlist(structure(x)$trans)
    glist = split(ulstruct, unlist(gid_by_exon))
    glengths = sapply(width(reduce(glist)), sum)
    tlengths = sapply(width(structure(x)$trans), sum)
    tfrags = lapply(1:nrow(tmeas), function(i){
        (tlengths[i]/1000) * tmeas[i,]
    }) ## still a bit slow
    tfrags = matrix(unlist(tfrags, use.names=FALSE), nrow=length(tfrags), byrow=TRUE)
    expr = t(sapply(1:length(inds_by_gene), function(i){colSums(tfrags[inds_by_gene[[i]],,drop=FALSE]) / glengths[i]}))
    rownames(expr) = names(inds_by_gene)
    return(expr)
})



