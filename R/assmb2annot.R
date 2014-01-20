#' match assembled transcript(s) to closest annotated transcript/gene
#' 
#' @param gtf path to a GTF file containing locations of annotated transcripts
#' @param chr which chromosome is your assembled transcript of interest on? Should be formatted the same way as chromosomes in \code{gtf}.
#' @param assembled GRanges object, with ranges representing the assembled transcript's exons. 
#' @return a GRangesList containing the annotated transcript(s) that most overlap \code{assembled}.
#' @author Alyssa Frazee
#' @export

assmb2annot = function(gtf, chr, assembled){
    # read in annotation and convert to GRangesList
    annot = gffRead(gtf)
    annotChr = subset(annot, feature=="exon" & seqname==chr)
    annotChr$transcript_id = getAttributeField(annotChr$attributes, "transcript_id")
    annotChr$transcript_id = sapply(annotChr$transcript_id, function(x) substr(x, 2, nchar(x)-1))
    annotChr$gene_id = getAttributeField(annotChr$attributes, "gene_id")
    annotChr$gene_id = sapply(annotChr$gene_id, function(x) substr(x, 2, nchar(x)-1))
    annotgr = split(GRanges(seqnames=Rle(chr),
        ranges=IRanges(start=annotChr$start, end=annotChr$end),
        strand=annotChr$strand, gene_id=annotChr$gene_id), annotChr$transcript_id)

    # make assembled transcript into GRangesList and make sure there's only one
    assembled = GRangesList(assembled)
    stopifnot(length(assembled) == 1)

    # find overlapping annotated transcripts:
    ol = findOverlaps(assembled, annotgr)

    # return nothing if there isn't any overlap:
    if(length(ol) == 0){
        message('Novel assembled transcript! (your assembled transcript does not overlap any annotated transcripts)')
        return(NULL)
    }

    # return one hit if only overlaps one:
    if(length(ol) == 1){
        return(annotgr[subjectHits(ol)])
    }

    # otherwise, find transcript with max % overlap and return that
    pctol = sapply(1:length(ol), function(i){
        # you need this so that the chromosomes for both GRanges match up :/ 
        t1 = annotgr[[subjectHits(ol)[i]]]
        t2 = assembled[[queryHits(ol)[i]]]
        t1chr = as.character(runValue(seqnames(t1)))
        t2chr = as.character(runValue(seqnames(t2)))
        t1good = GRanges(seqnames=Rle(t1chr), ranges=ranges(t1), strand=strand(t1))
        t2good = GRanges(seqnames=Rle(t2chr), ranges=ranges(t2), strand=strand(t2))
        return(pctOverlap(t1good, t2good))
    }) 
    return(annotgr[subjectHits(ol)[which(pctol == max(pctol))]])
}





