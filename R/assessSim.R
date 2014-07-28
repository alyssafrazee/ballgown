assessSim = function(bg, bgresults, annotation, chr, trulyDEids, 
    cuffdiffFile=NULL, qcut=0.05, UCSC=TRUE, ret=FALSE, nClosest=1, 
    limmaresults=NULL, plot=TRUE){
    
    stopifnot(all(qcut >= 0 & qcut <= 1))

    assemblygr = structure(bg)$trans
    if(class(annotation)!='data.frame'){
        annot = gffRead(annotation)
        annotsub = annot[annot[,3]=='exon' & annot[,1]==chr,] 
        #^ col3=feature, col1=seqname
        annotsub$tx = getAttributeField(annotsub$attributes, "transcript_id")
        if(UCSC){
            # strip quotes and strip off any "_2" business
            annotsub$tx = substr(annotsub$tx, 2, nchar(annotsub$tx)-1)
            annotsub$tx = sapply(annotsub$tx, 
                function(x){
                    paste(strsplit(x, split="_")[[1]][1:2],collapse="_")
                })
        }
    }else{
        annotsub = annotation
    }
    degtf = annotsub[annotsub$tx %in% trulyDEids,]
    degr = split(GRanges(seqnames=Rle(chr), 
        ranges=IRanges(start=degtf$start, end=degtf$end), 
        strand=degtf$strand), degtf$tx)

    deoverlaps = annotate_assembly(assemblygr, degr)
    ol_list = split(deoverlaps[, c(1, 3)], deoverlaps[,2])
    ## ties are handled according to "sort" -- 
    ## i.e., if there's a max overlap tie, the transcript 
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
    non_de_ids = setdiff(unique(texpr(bg, "all")$t_id), 
        as.numeric(unlist(truly_de_ids)))

    # load in cuffdiff results
    if(!is.null(cuffdiffFile)){
        cuff = read.table(cuffdiffFile, sep = "\t", header = TRUE)
        cuffok = cuff[cuff$status == 'OK',]
    }

    if(length(qcut) == 1){
        if(!is.null(cuffdiffFile)){
            cuff_decalls = as.character(cuffok$test_id[cuffok$q_value < qcut])
            cuff_decalls = texpr(bg, "all")$t_id[match(cuff_decalls, 
                texpr(bg, "all")$t_name)]
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
            cufffdr = sum(
                !(cuff_decalls %in% unlist(truly_de_ids)))/length(cuff_decalls)
        } else {
            cuffsens = cuffspec = cufffdr = NULL
        }
        
        if(!is.null(limmaresults)){ 
            limmasens = sum(found_by_limma) / length(found_by_limma)
            limmaspec = sum(
                !(non_de_ids %in% limma_decalls)) / length(non_de_ids)
            limmafdr = sum(!(limma_decalls %in% 
                unlist(truly_de_ids)))/length(limma_decalls)
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
            cuff_decalls = as.character(cuffok$test_id[
                cuffok$q_value < qcut[i]])
            cuff_decalls = texpr(bg,'all')$t_id[match(cuff_decalls, 
                texpr(bg,'all')$t_name)]
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
        bgfdr[i] = sum(
            !(bg_decalls %in% unlist(truly_de_ids)))/length(bg_decalls)

        if(!is.null(cuffdiffFile)){
            cuffsens[i] = sum(found_by_cuffdiff) / length(found_by_cuffdiff)    
            cuffspec[i] = sum(
                !(non_de_ids %in% cuff_decalls))/length(non_de_ids) 
            cufffdr[i] = sum(
                !(cuff_decalls %in% unlist(truly_de_ids)))/length(cuff_decalls)
        }
        
        if(!is.null(limmaresults)){
            limmasens[i] = sum(found_by_limma) / length(found_by_limma)
            limmaspec[i] = sum(
                !(non_de_ids %in% limma_decalls)) / length(non_de_ids)
            limmafdr[i] = sum(!(limma_decalls %in% 
                unlist(truly_de_ids)))/length(limma_decalls)
        }
        
    }

    if(plot){
        plot(1-bgspec, bgsens, col="dodgerblue", type="l", 
            xlab="false positive rate", ylab="true positive rate", lwd=2, 
            ylim=c(0,1))
        if(!is.null(cuffdiffFile)){
            lines(1-cuffspec, cuffsens, col="orange", lwd=2)
            legend('bottomright', lty=c(1,1), lwd=c(2,2), 
                col=c("dodgerblue", "orange"), c("ballgown", "cuffdiff"))
        }
    }

    if(ret){
        isDE = bgresults$id %in% unlist(truly_de_ids)
        return(list(ballgownsens=bgsens, cuffdiffsens=cuffsens, 
            limmasens=limmasens, ballgownspec=bgspec, cuffdiffspec=cuffspec, 
            limmaspec=limmaspec, bgfdr=bgfdr, cufffdr=cufffdr, 
            limmafdr=limmafdr, isDE=isDE))
    }
}
