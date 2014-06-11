#' plot annotated and assembled transcripts together
#'
#' @param tx.assembled a GRangesList object where the GRanges objects in the list represent sets of
#' exons comprising assembled transcripts
#' @param tx.annotated a GRangesList object where the GRanges objects in the list represent sets of 
#' exons comprising annotated transcripts
#' @param ind integer; index of \code{tx.annotated} specifying which annotated transcript to plot.  
#' All transcripts (assembled and annotated) overlapping \code{tx.annotated[[ind]]} will be plotted. 
#' Default 1.
#' @param main optional character string giving the title for the resulting plot.  Default: 
#' "Assembled and Annotated Transcripts"
#' 
#' @return No return value, but a plot is produced with annotated transcripts on the bottom panel 
#' (shaded in gray) and assembled transcripts on the top panel (shaded with diagonal lines).
#' 
#' @author Alyssa Frazee
#' 
#' @export
checkAssembledTx = function(tx.assembled, tx.annotated, ind=1, 
    main='Assembled and Annotated Transcripts'){
  
    ol = findOverlaps(tx.annotated, tx.assembled)
    ol.self = findOverlaps(tx.annotated, tx.annotated) #plot any overlapping transcripts also.
    inds.asmb = subjectHits(ol)[which(queryHits(ol)==ind)]
  
    if(length(inds.asmb)==0){
        warning("This annotated transcript was not covered by any assembled transcripts. 
            No plot will be produced", call.=FALSE)
        return(0)
    }
  
    # plot setup:
    par(mar=c(5,3,4,2))
    annot.list = lapply(subjectHits(ol.self)[which(queryHits(ol.self)==ind)], 
        function(i) tx.annotated[[i]])
    names(annot.list) = names(tx.annotated)[subjectHits(ol.self)[which(queryHits(ol.self)==ind)]]
    annot.df = as.data.frame(GRangesList(annot.list))
    asmbl.list = lapply(inds.asmb, function(i) tx.assembled[[i]])
    names(asmbl.list) = names(tx.assembled)[inds.asmb]
    asmbl.df = as.data.frame(GRangesList(asmbl.list))
    xax = seq(min(min(asmbl.df$start), min(annot.df$start)), max(max(asmbl.df$end),
        max(annot.df$end)), by=1)
    numtx = length(unique(annot.df$element))+length(unique(asmbl.df$element))
  
    # plot base:
    plot(xax, rep(0, length(xax)), ylim=c(0.5, numtx+1), type="n", xlab="genomic position", 
        yaxt="n", ylab="")
  
    # plot the annotated transcripts:
    txind = 0
    for(tx in unique(annot.df$element)){
        txind = txind+1 #counter
        gtsub = annot.df[annot.df$element==tx,]
        gtsub = gtsub[order(gtsub$start),]
        for(exind in 1:dim(gtsub)[1]){
            polygon(x=c(gtsub$start[exind], gtsub$start[exind], gtsub$end[exind], gtsub$end[exind]), 
                y=c(txind-0.4, txind+0.4, txind+0.4, txind-0.4), col="gray30")
            if(exind!=dim(gtsub)[1]){
                lines(c(gtsub$end[exind],gtsub$start[exind+1]),c(txind, txind), lty=2, col="gray50")
            }
        }
        st.compare = gtsub$start[-1]
        en.compare = gtsub$end[-length(gtsub$end)]
        diffs = st.compare-en.compare
        if(any(diffs<0)) warning(paste("overlapping exons in annotated transcript",tx), call.=FALSE)
    }
    txind = txind+0.5
    abline(h=txind+0.25)

    # plot the assembled transcripts
    for(tx in unique(asmbl.df$element)){
        txind = txind+1
        gtsub = asmbl.df[asmbl.df$element==tx,]
        gtsub = gtsub[order(gtsub$start),]
        for(exind in 1:dim(gtsub)[1]){
            polygon(x=c(gtsub$start[exind], gtsub$start[exind], gtsub$end[exind], gtsub$end[exind]), 
                y=c(txind-0.4, txind+0.4, txind+0.4, txind-0.4), col="black", density=15, angle=45)
            if(exind!=dim(gtsub)[1]){
                lines(c(gtsub$end[exind],gtsub$start[exind+1]),c(txind, txind), lty=2, col="gray50")
            }
        }
        st.compare = gtsub$start[-1]
        en.compare = gtsub$end[-length(gtsub$end)]
        diffs = st.compare-en.compare
        if(any(diffs<0)) warning(paste("overlapping exons in assembled transcript",tx), call.=FALSE)
    }
    title(main)
    axis(side=2, at=c(median(1:length(unique(annot.df$element))), 
        median(1:length(unique(asmbl.df$element)))+length(unique(annot.df$element))+0.5), 
        labels=c("annotated", "assembled"), tick=FALSE)
}
