#' analyze results of a two-group simulation experiment
#'
#' [documentation in progress]
#' 
#' @param bg ballgown object from the simulated data
#' @param bgresults data frame resulting from a call to \code{stattest(bg,...)}. Should be a transcript-level test (i.e., \code{feature="transcript"} in \code{stattest}).
#' @param annotation either a path to an annotation gtf file, or a data.frame of an annotation file that was already read in and subset to exons on the right chromosome.
#' @param chr 
#' @param trulyDEids 
#' @param cuffdiffFile 
#' @param qcut either a number between 0 and 1 to be used as the q-value significance cutoff, or a vector like \code{seq(0,1,by=0.01)} (i.e., ranging from 0 to 1 in even increments). 
#' @param UCSC 
#' @details \code{trulyDEids} should be the transcripts that were set to be differentially expressed, identified the SAME WAY AS THEY ARE in \code{annotation}. This is super important!!! 
#' 
#' Also: if \code{qcut} is a vector, \code{assessSim} returns sensitivities/specificities and creates an ROC plot. If \code{qcut} is a single number, sensitivity and specificity for both methods are returned, using \code{qcut} as a q-value significance cutoff.
#' @return creates ROC curve and returns sensitivities/specificities used for said ROC curve, comparing cuffdiff to ballgown (as used to create \code{bg}). 
#' @author Alyssa Frazee
#' @export

assessSim = function(bg, bgresults, annotation, chr, trulyDEids, cuffdiffFile, qcut=0.05, UCSC=TRUE, ret=FALSE){
    require(ballgown)
    require(GenomicRanges)

    stopifnot(all(qcut >= 0 & qcut <= 1))

    assemblygr = structure(bg)$trans
    if(class(annotation)!='data.frame'){
        annot = gffRead(annotation)
        annotsub = subset(annot, feature=="exon" & seqname==chr)        
        annotsub$tx = getAttributeField(annotsub$attributes, "transcript_id")
        if(UCSC){
            # strip quotes and strip off any "_2" business
            annotsub$tx = substr(annotsub$tx, 2, nchar(annotsub$tx)-1)
            annotsub$tx = sapply(annotsub$tx, function(x) paste(strsplit(x, split="_")[[1]][1:2],collapse="_"))
        }
    }else{
        annotsub = annotation
    }
    
    degtf = subset(annotsub, tx %in% trulyDEids)
    degr = split(GRanges(seqnames = Rle(chr), ranges = IRanges(start = degtf$start, end = degtf$end), strand = degtf$strand), degtf$tx)

    deoverlaps = annotate_assembly(assemblygr, degr)
    ol_list = split(deoverlaps[, c(1, 3)], deoverlaps[,2])
    find_correct = function(x) {
        return(x[which(x[, 2] == max(x[, 2])), 1])
    }
    truly_de = lapply(ol_list, find_correct)

    # if there are ties in "closest assembled transcript", count the transcript correctly called DE if at least 1 of the closest ones is found
    inds_to_ids = function(x){
        as.numeric(names(assemblygr)[x])
    }##in case txid and indices do not match
    truly_de_ids = lapply(truly_de, inds_to_ids)
    non_de_ids = setdiff(unique(texpr(bg, "all")$t_id), as.numeric(unlist(truly_de_ids)))

    # load in cuffdiff results
    cuff = read.table(cuffdiffFile, sep = "\t", header = TRUE)
    cuffok = subset(cuff, status == "OK")
    
    if(length(qcut) == 1){
        cuff_decalls = as.character(cuffok$test_id[cuffok$q_value < qcut])
        cuff_decalls = texpr(bg, "all")$t_id[match(cuff_decalls, texpr(bg, "all")$t_name)]
        bg_decalls = bgresults$id[which(bgresults$qval < qcut)]
        bg_decalls = as.numeric(as.character(bg_decalls))
        found_by_bg = sapply(truly_de_ids, function(x){
            any(x %in% bg_decalls)
        })
        found_by_cuffdiff = sapply(truly_de_ids, function(x){
            any(x %in% cuff_decalls)
        })
        bgsens = sum(found_by_bg) / length(found_by_bg)
        cuffsens = sum(found_by_cuffdiff) / length(found_by_cuffdiff)
        bgspec = sum(!(non_de_ids %in% bg_decalls))/length(non_de_ids)
        cuffspec = sum(!(non_de_ids %in% cuff_decalls))/length(non_de_ids)
        bgfdr = sum(!(bg_decalls %in% unlist(truly_de_ids)))/length(bg_decalls)
        cufffdr = sum(!(cuff_decalls %in% unlist(truly_de_ids)))/length(cuff_decalls)
        return(list(ballgownsens = bgsens, cuffdiffsens = cuffsens, 
            ballgownspec = bgspec, cuffdiffspec = cuffspec, ballgownfdr = bgfdr, 
            cuffdifffdr = cufffdr))
    }

    # else, make ROC plots
    bgsens = cuffsens = bgspec = cuffspec = bgfdr = cufffdr = NULL
    for(i in 1:length(qcut)){
        cuff_decalls = as.character(cuffok$test_id[cuffok$q_value < qcut[i]])
        cuff_decalls = texpr(bg,'all')$t_id[match(cuff_decalls, texpr(bg,'all')$t_name)]
        bg_decalls = bgresults$id[which(bgresults$qval < qcut[i])]
        bg_decalls = as.numeric(as.character(bg_decalls))
        found_by_bg = sapply(truly_de_ids, function(x){
            any(x %in% bg_decalls)
        })
        found_by_cuffdiff = sapply(truly_de_ids, function(x){
            any(x %in% cuff_decalls)
        })
        bgsens[i] = sum(found_by_bg) / length(found_by_bg)
        cuffsens[i] = sum(found_by_cuffdiff) / length(found_by_cuffdiff)
        bgspec[i] = sum(!(non_de_ids %in% bg_decalls))/length(non_de_ids)
        cuffspec[i] = sum(!(non_de_ids %in% cuff_decalls))/length(non_de_ids) 
        bgfdr[i] = sum(!(bg_decalls %in% unlist(truly_de_ids)))/length(bg_decalls)
        cufffdr[i] = sum(!(cuff_decalls %in% unlist(truly_de_ids)))/length(cuff_decalls)
    }
    plot(1-bgspec, bgsens, col="dodgerblue", type="l", xlab="false positive rate", ylab="true positive rate", lwd=2, ylim=c(0,1))
    lines(1-cuffspec, cuffsens, col="orange", lwd=2)
    legend('bottomright', lty=c(1,1), lwd=c(2,2), col=c("dodgerblue", "orange"), c("ballgown", "cuffdiff"))
    if(ret){
        return(list(ballgownsens=bgsens, cuffdiffsens=cuffsens, 
            ballgownspec=bgspec, cuffdiffspec=cuffspec, bgfdr=bgfdr, cufffdr=cufffdr))
    }
}