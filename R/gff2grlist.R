#'Read in transcript file (gtf format) and return a GRangesList of transcripts
#'
#'Read a gtf file into memory, return a GRangesList of transcripts in that
#'file.  Each element (transcript) of the returned GRangesList is a GRanges object 
#'specifying the exons that comprise that transcript. Usually used to read in annotated
#'transcripts so that their format will match the format of assembled transcripts in 
#'a ballgown object (in \code{structure(ballgown.object)$trans}).
#'
#'
#'@param gtf string specifying location of gtf file to read
#'@param txidentifier string specifying the identifier in the of the "attributes" column 
#'of the the gtf file that gives the transcript id for each exon.
#'@param exonsonly if TRUE, return transcripts only as a set of exons - do not include CDS, start/stop codons, etc.
#'@return A \code{GRangesList} object with each GRangesElement denoting a transcript in the gtf file.
#'@author Alyssa Frazee
#'@export
gff2grlist = function(gtf, txidentifier = "transcript_id", exonsonly = TRUE){
  gtf.dataframe = gffRead(gtf)
  gtf.dataframe$txid = getAttributeField(gtf.dataframe$attributes, txidentifier, attrsep = "; ")
  gtf.dataframe = gtf.dataframe[,-9]
  if(exonsonly) gtf.dataframe = subset(gtf.dataframe, feature=="exon")
  gtf.list = split(gtf.dataframe, gtf.dataframe$txid)
  gtf.grl = lapply(gtf.list, function(x){
    x$strand[x$strand=="."] = "*"
    GRanges(seqnames=x$seqname, ranges=IRanges(start=x$start, end=x$end), strand = x$strand, id=rep(0, length(x$seqname)), transcripts=rep("NA", length(x$seqname)))
  })
  return(GRangesList(gtf.grl))
}