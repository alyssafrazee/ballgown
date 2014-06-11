#' analyze results of a two-group simulation experiment
#'
#' This function calculates sensitivity, specificity, and false discovery rates for a simulation 
#' experiment. Currently it compares ballgown, cuffdiff, and limma. Assembled transcripts in the 
#' ballgown objects are matched to annotated transcripts using percent overlap as distance. Analysis
#' is done one chromosome at a time.
#' 
#' @param bg ballgown object from the simulated data
#' @param bgresults data frame resulting from a call to \code{stattest(bg,...)}. Should be a 
#' transcript-level test (i.e., \code{feature="transcript"} in \code{stattest}).
#' @param annotation either a path to an annotation gtf file, or a data.frame of an annotation file 
#' that was already read in and subset to exons on the right chromosome.
#' @param chr Chromosome to analyze. Currently you can only do one at a time. 
#' @param trulyDEids character vector of transcript IDs (matching the IDs in \code{annotation}) that
#' were set to be truly differentially expressed in the simulation
#' @param cuffdiffFile path to the Cuffdiff output file \code{isoform_exp.diff}.
#' @param qcut either a number between 0 and 1 to be used as the q-value significance cutoff, or a 
#' vector like \code{seq(0,1,by=0.01)} (i.e., ranging from 0 to 1 in even increments). 
#' @param UCSC set to \code{TRUE} if \code{annotation} is UCSC annotation (UCSC annotations require)
#' extra text processing
#' @param ret set to TRUE to return simulation information in addition to plotting ROC curve. (FALSE
#' by default).
#' @param nClosest for true positives, make a positive call if any of the \code{nClosest} closest 
#' (most overlapping) transcripts are called DE by a method. Default 1.
#' @param limmaresults data frame or list with transcript q-values in one column ("qval") and 
#' transcript ids (matching those in \code{bg}) in another column ("id").
#' @param plot set to FALSE to disable the plotting of the ROC curve when \code{qcut} is a vector. 
#' (TRUE by default).
#' 
#' @details If \code{qcut} is a vector, \code{assessSim} returns sensitivities/specificities and 
#' creates an ROC plot, if \code{plot} is TRUE. If \code{qcut} is a single number, sensitivity and 
#' specificity for both methods are returned, using \code{qcut} as a q-value significance cutoff.
#' 
#' @return creates ROC curve and returns sensitivities/specificities used for said ROC curve, 
#' comparing ballgown (as used to create \code{bg}) to cuffdiff and limma.
#' 
#' @author Alyssa Frazee
#' 
#' @export

assessSim = function(bg, bgresults, annotation, chr, trulyDEids, 
    cuffdiffFile=NULL, qcut=0.05, UCSC=TRUE, ret=FALSE, nClosest=1, 
    limmaresults=NULL, plot=TRUE){
    
    stopifnot(all(qcut >= 0 & qcut <= 1))

    assemblygr = structure(bg)$trans
    if(class(annotation)!='data.frame'){
        annot = gffRead(annotation)
        annotsub = annot[annot[,3]=='exon' & annot[,1]==chr,] #col3=feature, col1=seqname
        annotsub$tx = getAttributeField(annotsub$attributes, "transcript_id")
        if(UCSC){
            # strip quotes and strip off any "_2" business
            annotsub$tx = substr(annotsub$tx, 2, nchar(annotsub$tx)-1)
            annotsub$tx = sapply(annotsub$tx, 
                function(x) paste(strsplit(x, split="_")[[1]][1:2],collapse="_"))
        }
    }else{
        annotsub = annotation
    }
    degtf = annotsub[annotsub$tx %in% trulyDEids,]
    degr = split(GRanges(seqnames=Rle(chr), 
        ranges=IRanges(start=degtf$start, end=degtf$end), strand=degtf$strand), degtf$tx)

    deoverlaps = annotate_assembly(assemblygr, degr)
    ol_list = split(deoverlaps[, c(1, 3)], deoverlaps[,2])
    ## ties are handled according to "sort" -- i.e., if there's a max overlap tie, the transcript 
    ## with the smaller start position will be chosen
    find_correct = function(x){
        maxInd = min(nrow(x), nClosest)
        return(x[order(x[,2], decreasing=TRUE)[1:maxInd] ,1])
    }
    truly_de = lapply(ol_list, find_correct)

    inds_to_ids = function(x){
        as.numeric(names(assemblygr)[x])
    }##in case txid and indices do not match
    truly_de_ids = lapply(truly_de, inds_to_ids)
    non_de_ids = setdiff(unique(texpr(bg, "all")$t_id), as.numeric(unlist(truly_de_ids)))

    # load in cuffdiff results
    if(!is.null(cuffdiffFile)){
        cuff = read.table(cuffdiffFile, se p = "\t", header = TRUE)
        cuffok = cuff[cuff$status == 'OK',]
    }

    if(length(qcut) == 1){
        if(!is.null(cuffdiffFile)){
            cuff_decalls = as.character(cuffok$test_id[cuffok$q_value < qcut])
            cuff_decalls = texpr(bg, "all")$t_id[match(cuff_decalls, texpr(bg, "all")$t_name)]
            found_by_cuffdiff = sapply(truly_de_ids, function(x){
                any(x %in% cuff_decalls)
            })
        }
        bg_decalls = bgresults$id[which(bgresults$qval < qcut)]
        bg_decalls = as.numeric(as.character(bg_decalls))
        found_by_bg = sapply(truly_de_ids, function(x){
            any(x %in% bg_decalls)
        })
        if(!is.null(limmaresults)) { 
            limma_decalls = limmaresults$id[limmaresults$qval < qcut]
            found_by_limma = sapply(truly_de_ids, function(x){
                any(x %in% limma_decalls)
            })
        }

        bgsens = sum(found_by_bg) / length(found_by_bg)
        bgspec = sum(!(non_de_ids %in% bg_decalls))/length(non_de_ids)
        bgfdr = sum(!(bg_decalls %in% unlist(truly_de_ids)))/length(bg_decalls)

        if(!is.null(cuffdiffFile)){
            cuffsens = sum(found_by_cuffdiff) / length(found_by_cuffdiff)    
            cuffspec = sum(!(non_de_ids %in% cuff_decalls))/length(non_de_ids)
            cufffdr = sum(!(cuff_decalls %in% unlist(truly_de_ids)))/length(cuff_decalls)
        } else {
            cuffsens = cuffspec = cufffdr = NULL
        }
        
        if(!is.null(limmaresults)){ 
            limmasens = sum(found_by_limma) / length(found_by_limma)
            limmaspec = sum(!(non_de_ids %in% limma_decalls)) / length(non_de_ids)
            limmafdr = sum(!(limma_decalls %in% unlist(truly_de_ids)))/length(limma_decalls)
        } else {
            limmasens = limmaspec = limmafdr = NULL
        }
        
        return(list(ballgownsens=bgsens, cuffdiffsens=cuffsens, 
            limmasens=limmasens, ballgownspec=bgspec, 
            cuffdiffspec=cuffspec, limmapsec=limmaspec, 
            ballgownfdr=bgfdr, cuffdifffdr=cufffdr, limmafdr=limmafdr))
    }

    # else, make ROC plots
    bgsens = cuffsens = bgspec = cuffspec = bgfdr = cufffdr = NULL
    limmasens = limmaspec = limmafdr = NULL
    for(i in 1:length(qcut)){

        if(!is.null(cuffdiffFile)){
            cuff_decalls = as.character(cuffok$test_id[cuffok$q_value < qcut[i]])
            cuff_decalls = texpr(bg,'all')$t_id[match(cuff_decalls, texpr(bg,'all')$t_name)]
            found_by_cuffdiff = sapply(truly_de_ids, function(x){
                any(x %in% cuff_decalls)
            })
        }

        bg_decalls = bgresults$id[which(bgresults$qval < qcut[i])]
        bg_decalls = as.numeric(as.character(bg_decalls))
        found_by_bg = sapply(truly_de_ids, function(x){
            any(x %in% bg_decalls)
        })

        if(!is.null(limmaresults)) { 
            limma_decalls = limmaresults$id[limmaresults$qval < qcut[i]]
            found_by_limma = sapply(truly_de_ids, function(x){
                any(x %in% limma_decalls)
            })
        }

        bgsens[i] = sum(found_by_bg) / length(found_by_bg)
        bgspec[i] = sum(!(non_de_ids %in% bg_decalls))/length(non_de_ids)
        bgfdr[i] = sum(!(bg_decalls %in% unlist(truly_de_ids)))/length(bg_decalls)

        if(!is.null(cuffdiffFile)){
            cuffsens[i] = sum(found_by_cuffdiff) / length(found_by_cuffdiff)    
            cuffspec[i] = sum(!(non_de_ids %in% cuff_decalls))/length(non_de_ids) 
            cufffdr[i] = sum(!(cuff_decalls %in% unlist(truly_de_ids)))/length(cuff_decalls)
        }
        
        if(!is.null(limmaresults)){
            limmasens[i] = sum(found_by_limma) / length(found_by_limma)
            limmaspec[i] = sum(!(non_de_ids %in% limma_decalls)) / length(non_de_ids)
            limmafdr[i] = sum(!(limma_decalls %in% unlist(truly_de_ids)))/length(limma_decalls)
        }
        
    }

    if(plot){
        plot(1-bgspec, bgsens, col="dodgerblue", type="l", xlab="false positive rate", 
            ylab="true positive rate", lwd=2, ylim=c(0,1))
        if(!is.null(cuffdiffFile)){
            lines(1-cuffspec, cuffsens, col="orange", lwd=2)
            legend('bottomright', lty=c(1,1), lwd=c(2,2), col=c("dodgerblue", "orange"), 
                c("ballgown", "cuffdiff"))
        }
    }

    if(ret){
        isDE = bgresults$id %in% unlist(truly_de_ids)
        return(list(ballgownsens=bgsens, cuffdiffsens=cuffsens, limmasens=limmasens,
            ballgownspec=bgspec, cuffdiffspec=cuffspec, limmaspec=limmaspec,
            bgfdr=bgfdr, cufffdr=cufffdr, limmafdr=limmafdr, isDE=isDE))
    }
}
