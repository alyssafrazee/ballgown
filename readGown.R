##### READGOWN + HELPERS
##### CREDIT FOR ORIGINAL goes to Leonardo - we can merge into one "publishable" repository later
##### updates made to his readGown function:
##### changed some of the messages
##### changed the "." strand identifier to the "*" strand identifier.
##### added GRanges & pData (+ some additional) objects to the output

#### Main function
readGown <- function(dataDir, samplePattern, bamfiles = NULL, pData = NULL, verbose=TRUE, ...) {
  if(verbose) message(date())

  ## Load required pkgs
  suppressMessages(library(plyr))
  suppressMessages(library(GenomicRanges))

  if(FALSE){
    # for testing
    dataDir <- "/amber2/scratch/lcollado/brain_rna/ballgown"
    samplePattern <- "orb"
  }

  ## Identify the sample directories
  dirs <- list.files(path=dataDir, pattern=samplePattern, full.names=TRUE)
  names(dirs) <- list.files(path=dataDir, pattern=samplePattern)
  n <- length(dirs)

  if(FALSE){
    # for testing
    n <- 3
  }

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
   t2g = data.frame(t_id = trans$t_id, g_id = as.numeric(sapply(trans$gene_id, function(s) substr(s, 6, nchar(s)))))

   ## Read phenotype table, if given:
   if(!is.null(pData)){
       if(verbose) message(paste0(date(),": Reading phenotype table"))
       phx = read.table(pData, stringsAsFactors=FALSE, ...)
       theorder = sapply(names(dirs), function(x) which(phx$dirname==x))
       phx = phx[theorder,]
   }
   if(is.null(pData)) phx = NULL

  if(verbose) message("Wrapping up the results")
  result <- list(data = list(intron=intron, exon=exon, trans=trans), indexes=list(e2t=e2t, i2t=i2t, t2g=t2g, bamfiles = bamfiles, pData = phx), structure = list(intron = introngr, exon = exongr, trans = transgrl), dirs=dirs, mergedDate=date())

  if(verbose) message(date())
  # Done!
  return(result)
}

#### Auxiliary functions

## Read intron files
.readIntron <- function(file){
  intron <- read.table(file, header=TRUE, sep="\t", colClasses=c("integer", "character", "factor", "integer", "integer", "integer", "integer", "numeric"))
  intron <- intron[order(intron$i_id), ]
  rownames(intron) <- 1:nrow(intron)
  return(intron)
}

## Read counts and raw coverage
.readExon <- function(file) {
  exon <- read.table(file, header=TRUE, sep="\t", colClasses=c("integer", "character", "factor", "integer", "integer", "integer", "integer", "numeric", "numeric", "numeric", "numeric", "numeric"))
  exon <- exon[order(exon$e_id), ]
  rownames(exon) <- 1:nrow(exon)
  return(exon)
}

## Read transcript data files
.readTrans <- function(file) {
  trans <- read.table(file, header=TRUE, sep="\t", colClasses=c("integer", "character", "factor", "integer", "integer", "character", "integer", "integer", "character", "character", "numeric", "numeric"))
  trans <- trans[order(trans$t_id), ]
  rownames(trans) <- 1:nrow(trans)
  return(trans)
}


testobj = readGown(dataDir = "/amber2/scratch/lcollado/brain_rna/ballgown", samplePattern = "orb", sep="\t", header = TRUE)
save(testobj, file="/amber2/scratch/jleek/orbFrontal/results/ballgown_structure_test.rda")

