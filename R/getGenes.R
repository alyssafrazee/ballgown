#' label assembled transcripts with gene names
#' 
#' @param gtf path to a GTF file containing locations of annotated transcripts
#' @param assembled GRangesList object, with each set of ranges representing exons of an assembled transcript.
#' @param UCSC set to \code{TRUE} if you're using a UCSC gtf file. (Requires some extra text-processing).
#' @param attribute set to gene_id (default) if you want the gene ID or gene_name if you want the gene symbol 
#' @return a character vector of the same length as \code{assembled}, providing the name(s) of the gene(s) that overlaps each transcript in \code{assembled}.
#' @details chromosome labels in \code{gtf} and \code{assembled} should match. (i.e., you should provide the path to a gtf corrsponding to the same annotation you used when constructing \code{assembled})
#' 
#' @author Alyssa Frazee
#' @export

getGenes = function(gtf, assembled, UCSC=TRUE, attribute = "gene_id"){
    # read in annotation and split by gene:
	if(!attribute %in%  c("gene_id", "gene_name")) stop("attribute must be gene_id or gene_name\n")
    require(GenomicRanges)
    annot = gffReadGR(gtf)
    annotEx = annot[mcols(annot)$type=="exon"]
    annotEx$gene_id = getAttributeField(as.character(mcols(annotEx)$group), attribute)
    if(UCSC){
        # strip quotes off of gene names
        annotEx$gene_id = substr(annotEx$gene_id, 2, nchar(annotEx$gene_id)-1)
    }
    geneLocs = split(annotEx, annotEx$gene_id)

    # find which transcripts overlap which genes:
    ol = findOverlaps(assembled, geneLocs)
    split_ol = split(names(geneLocs[subjectHits(ol)]), queryHits(ol)) #to deal with transcripts overlapping >1 gene
	split_ol = sapply(split_ol, function(x) paste(x, collapse="; "))
    gene_id = rep("", length(assembled))
    gene_id[as.numeric(names(split_ol))] = as.character(split_ol)

	gene_id = CharacterList(strsplit(gene_id,"; "))
		
    return(gene_id)
}


