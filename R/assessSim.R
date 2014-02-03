#' analyze results of a two-group simulation experiment
#'
#' @param bg ballgown object from the simulated data
#' @param bgresults data frame resulting from a call to \code{stattest(bg,...)}. Should be a transcript-level test (i.e., \code{feature="transcript"} in \code{stattest}).
#' @param annotation
#' @param chr 
#' @param trulyDEids 
#' @param cuffdiffFile 
#' @param qcut either a number between 0 and 1 to be used as the q-value significance cutoff, or a vector like \code{seq(0,1,by=0.01)} (i.e., ranging from 0 to 1 in even increments). 
#' @param UCSC 
#' @details 
#' @return 
#' @author Alyssa Frazee
#' @export

assessSim = function(bg, bgresults, annotation, chr, trulyDEids, cuffdiffFile, qcut=0.05, UCSC=TRUE){
    require(ballgown)
    require(GenomicRanges)

    #trulyDEids should match transcript IDs in the "annotation" file.
    # if qcut is a vector, this will return things that are useful for ROC plots

    assemblygr = structure(bg)$trans
    annot = gffRead(annotation)
    annotsub = subset(annot, feature=="exon" & seqname==chr)        
    annotsub$tx = getAttributeField(annotsub$attributes, "transcript_id")
    if(UCSC){
        # strip quotes and strip off any "_2" business
        annotsub$tx = sapply(annotsub$tx, function(x) substr(x,2,nchar(x)-1))
        annotsub$tx = sapply(annotsub$tx, function(x) paste(strsplit(x, split="_")[[1]][1:2],collapse="_"))
    }
    
    degtf = subset(annotsub, tx %in% trulyDEids)
    stopifnot(length(unique(degtf$tx)) == length(unique(trulyDEids)))
    nondegtf = subset(annotsub, !(tx %in% trulyDEids))

    degr = split(GRanges(seqnames=Rle(chr),
        ranges=IRanges(start=degtf$start, end=degtf$end),
        strand=degtf$strand), degtf$tx)
    nondegr = split(GRanges(seqnames=Rle(chr),
        ranges=IRanges(start=nondegtf$start, end=nondegtf$end),
        strand=nondegtf$strand), nondegtf$tx)

    deOverlaps = findOverlaps(assemblygr, degr)
    nondeOverlaps = findOverlaps(assemblygr, nondegr)

    depercent = function(i){
        t1 = degr[[subjectHits(deOverlaps)[i]]]
        t2 = assemblygr[[queryHits(deOverlaps)[i]]]
        t1chr = as.character(runValue(seqnames(t1)))
        t2chr = as.character(runValue(seqnames(t2)))
        t1good = GRanges(seqnames=Rle(t1chr), ranges=ranges(t1), strand=strand(t1))
        t2good = GRanges(seqnames=Rle(t2chr), ranges=ranges(t2), strand=strand(t2))
        return(pctOverlap(t1good, t2good))
    }

    nondepercent = function(i){
        t1 = nondegr[[subjectHits(nondeOverlaps)[i]]]
        t2 = assemblygr[[queryHits(nondeOverlaps)[i]]]
        t1chr = as.character(runValue(seqnames(t1)))
        t2chr = as.character(runValue(seqnames(t2)))
        t1good = GRanges(seqnames=Rle(t1chr), ranges=ranges(t1), strand=strand(t1))
        t2good = GRanges(seqnames=Rle(t2chr), ranges=ranges(t2), strand=strand(t2))
        return(pctOverlap(t1good, t2good))
    }

    pctDE = sapply(1:length(deOverlaps), depercent)
    #pctNonDE = sapply(1:length(nondeOverlaps), nondepercent)


    corresponding_to_de = data.frame(
        assembled=queryHits(deOverlaps), 
        annotated=subjectHits(deOverlaps), 
        overlap=pctDE)
    if(any(!(1:length(trulyDEids) %in% corresponding_to_de$annotated))){
        message('some DE transcripts did not overlap any assembled transcripts. Indexes:')
        print(which(!(1:length(trulyDEids) %in% corresponding_to_de$annotated)))
    }

    ol_list = split(corresponding_to_de[,c(1,3)], 
        corresponding_to_de$annotated)
    find_correct = function(x){
        return(x[which(x[,2]==max(x[,2])),1])
    }
    truly_de = lapply(ol_list, find_correct)
    num_max = sapply(truly_de, length)
    stopifnot(all(num_max == 1))

    bgtxids = texpr(bg,'all')$t_name[match(bgresults$id, texpr(bg,'all')$t_id)]

    ### BALLGOWN IDS that should and should not be DE:
    de_ids = unique(as.numeric(names(assemblygr)[as.numeric(truly_de)]))
    non_de_ids = setdiff(unique(texpr(bg,'all')$t_id), as.numeric(de_ids))


    # load in cuffdiff results:
    cuff = read.table(cuffdiffFile, sep='\t', header=TRUE)
    cuffok = subset(cuff, status=='OK') 

    if(length(qcut) == 1){
        cuff_decalls = as.character(cuffok$test_id[cuffok$q_value < qcut])
        # translate cufflinks DE calls into ballgown transcript ids (numeric):
        cuff_decalls = texpr(bg,'all')$t_id[match(cuff_decalls, texpr(bg,'all')$t_name)]
        bg_decalls = bgresults$id[which(bgresults$qval < qcut)]
        bg_decalls = as.numeric(as.character(bg_decalls))

        ### sensitivity
        bgsens = sum(de_ids %in% bg_decalls)/length(de_ids)
        cuffsens = sum(de_ids %in% cuff_decalls)/length(de_ids)

        ### specificity
        bgspec = sum(!(non_de_ids %in% bg_decalls))/length(non_de_ids)
        cuffspec = sum(!(non_de_ids %in% cuff_decalls))/length(non_de_ids) 

        ### false discovery rates
        bgfdr = sum(!(bg_decalls %in% de_ids))/length(bg_decalls)
        cufffdr = sum(!(cuff_decalls %in% de_ids))/length(cuff_decalls) 

        return(list(ballgownsens=bgsens, cuffdiffsens=cuffsens, 
            ballgownspec=bgspec, cuffdiffspec=cuffspec,
            ballgownfdr=bgfdr, cuffdifffdr=cufffdr))
        }
    ## else, make ROC plots:
    bgsens = cuffsens = bgspec = cuffspec = NULL
    for(i in 1:length(qcut)){
        cuff_decalls = as.character(cuffok$test_id[cuffok$q_value < qcut[i]])
        cuff_decalls = texpr(bg,'all')$t_id[match(cuff_decalls, texpr(bg,'all')$t_name)]
        bg_decalls = bgresults$id[which(bgresults$qval < qcut[i])]
        bgsens[i] = sum(de_ids %in% bg_decalls)/length(de_ids)
        cuffsens[i] = sum(de_ids %in% cuff_decalls)/length(de_ids)
        bgspec[i] = sum(!(non_de_ids %in% bg_decalls))/length(non_de_ids)
        cuffspec[i] = sum(!(non_de_ids %in% cuff_decalls))/length(non_de_ids) 
    }
    plot(1-bgspec, bgsens, col="dodgerblue", type="l", xlab="false positive rate", ylab="true positive rate", lwd=2, ylim=c(0,1))
    lines(1-cuffspec, cuffsens, col="orange", lwd=2)
    return(list(ballgownsens=bgsens, cuffdiffsens=cuffsens, 
        ballgownspec=bgspec, cuffdiffspec=cuffspec))

}

# cuffFile = "~/Desktop/simcuff.txt"
# bgresults = stattest(simgown, feature="transcript", meas="FPKM", covariate="group")
# annFile = "~/Google Drive/hopkins/research/_cufflinks visualization project/genes-clean-small.gtf"
# de1 = gffRead('~/Google Drive/hopkins/research/_cufflinks visualization project/de-tx1.gtf')
# de1 = subset(de1, feature=="exon")
# de2 = gffRead('~/Google Drive/hopkins/research/_cufflinks visualization project/de-tx2.gtf')
# de2 = subset(de2, feature=="exon")
# de1$tx = getAttributeField(de1$attributes, "transcript_id")
# de2$tx = getAttributeField(de2$attributes, "transcript_id")
# trulyDEids = unique(c(de1$tx,de2$tx))

# assessSim(bg=simgown, 
#     bgresults=bgresults, 
#     annotation=annFile, 
#     chr="22", 
#     trulyDEids=trulyDEids,
#     cuffdiffFile=cuffFile)


#assessSim(bg=simgown, bgresults=bgresults, annotation=annFile, 
 #    chr="22", trulyDEids=trulyDEids, cuffdiffFile=cuffFile,UCSC=TRUE)
# assessSim(bg=simgown, bgresults=bgresults, annotation=annFile, 
#      chr="22", trulyDEids=trulyDEids, cuffdiffFile=cuffFile, 
#      qcut=seq(0,1,by=0.01),UCSC=FALSE)


