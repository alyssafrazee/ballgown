#' match assembled transcript(s) to closest annotated transcript/gene
#' 
#' @param gtf path to a GTF file containing locations of annotated transcripts
#' @param chr which chromosome is your assembled transcript of interest on? Should be formatted the same way as chromosomes in \code{gtf}.
#' @param assembled GRangesList object, with each set of ranges representing exons of an assembled transcript.
#' @return a GRangesList containing the corresponding annotated transcript(s) that most overlap \code{assembled}.
#' @details the \code{elementMetadata} slot of each \code{GRanges} object in the returned \code{GRangesList} contains the annotated \code{gene_id} and \code{transcript_id}. The \code{names} component of this \code{GRangesList} gives the assembled transcript ids. 
#' 
#' Also be careful not to confuse this with \code{annot2assmb}, which finds the closest \emph{assembled} transcript to each \emph{annotated} transcript. That function is more useful in simulations; this one is more useful for getting biological information out of an assembly.
#' @author Alyssa Frazee
#' @export

assmb2annot = function(gtf, assembled){
    # read in annotation and convert to GRangesList
    require(GenomicRanges)
    annot = gffRead(gtf)
    annotEx = subset(annot, feature=="exon")
    annotEx$transcript_id = getAttributeField(annotEx$attributes, "transcript_id")
    annotEx$transcript_id = substr(annotEx$transcript_id, 2, nchar(annotEx$transcript_id)-1))
    annotEx$gene_id = getAttributeField(annotEx$attributes, "gene_id")
    annotEx$gene_id = substr(annotEx$gene_id, 2, nchar(annotEx$gene_id)-1))
    annotgr = split(GRanges(seqnames=Rle(annotEx$seqname),
        ranges=IRanges(start=annotEx$start, end=annotEx$end),
        strand=annotEx$strand, gene_id=annotEx$gene_id), annotEx$transcript_id)
    message('done reading gtf file')

    # find overlapping annotated transcripts:
    ol = findOverlaps(assembled, annotgr)

    # return nothing if there isn't any overlap:
    if(length(ol) == 0){
        message('All assembled transcripts are novel! That\'s fun.')
        return(NULL)
    }

    # return one hit if only overlaps one:
    if(length(ol) == 1){
        ret = annotgr[subjectHits(ol)]
        elementMetadata(ret[[1]])$asmbtx_id = rep(names(assembled)[queryHits(ol)], nrow(elementMetadata(ret[[1]])))
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
    allInfo = data.frame(assembled=queryHits(ol), annotated=subjectHits(ol), overlap=pctol)
    olList = split(allInfo[,c(2,3)], allInfo$assembled)
    olTx = lapply(olList, function(x) x[which(x[,2] == max(x[,2])), 1])
    ret = annotgr[as.numeric(olTx)]
    ann_names = names(ret)
    ret = lapply(seq_along(ret), function(i){
        ambInd = as.numeric(names(olTx)[i])
        addedEMD = ret[[i]]
        elementMetadata(addedEMD)$transcript_id = rep(ann_names[i], length(ret[[i]]))
        return(addedEMD)
    })
    names(ret) = names(assembled)[as.numeric(names(olList))]
    return(GRangesList(ret))
}




