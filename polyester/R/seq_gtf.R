#' @title Get transcript sequences from GTF file and sequence info
#'
#' @description Given a GTF file (for transcript structure) and DNA sequences, return a DNAStringSet
#' of transcript sequences
#' @param gtf path to GTF file
#' @param seqpath path to folder containing one FASTA file (\code{.fa} extension) for each 
#' chromosome in \code{gtf}
#' @param exononly if \code{TRUE} (as it is by default), only create transcript sequences from the 
#' features labeled \code{exon} in \code{gtf}
#' @param idfield in the \code{attributes} column of \code{gtf}, what is the name of the field 
#' identifying transcripts? Should be character. Default \code{"transcript_id"}.
#' @param attrsep in the \code{attributes} column of \code{gtf}, how are attributes separated? 
#' Default \code{"; "}.
#' @export
#' @references \url{http://www.ensembl.org/info/website/upload/gff.html}
#' @return DNAStringSet containing transcript sequences, with names corresponding to \code{idfield}
#' in \code{gtf}
#' @importFrom ballgown gffRead
#' @importFrom ballgown getAttributeField
seq_gtf = function(gtf, seqpath, exononly=TRUE, idfield="transcript_id", attrsep="; "){
    gtf_dat = gffRead(gtf)
    if(exononly){
        gtf_dat = subset(gtf_dat, feature=="exon")
    }

    chrs = unique(gtf_dat$seqname)
    fafiles = list.files(seqpath)
    if(!(all(paste0(chrs, '.fa') %in% fafiles))){
        stop("all chromosomes in the GTF file must have .fa files in seqpath")
    }

    seqlist = lapply(chrs, function(chr){
        dftmp = subset(gtf_dat, seqname==chr)
        fullseq = readDNAStringSet(paste0(seqpath, '/', chr, '.fa'))
        these_seqs = subseq(rep(fullseq, times=nrow(dftmp)), start=dftmp$start, end=dftmp$end)
        names(these_seqs) = getAttributeField(dftmp$attributes, idfield, attrsep=attrsep)
        these_seqs
    })

    full_list = do.call(c, seqlist)
    split_list = split(full_list, names(full_list))
    DNAStringSet(lapply(split_list, unlist)) #took 340 sec on whole human transcriptome hg19
}

