#' plot annotated and assembled transcripts together
#'
#' @param assembled a GRangesList object where the GRanges objects in the list 
#'   represent sets of exons comprising assembled transcripts
#' @param annotated a GRangesList object where the GRanges objects in the list
#'   represent sets of 
#' exons comprising annotated transcripts
#' @param ind integer; index of \code{annotated} specifying which annotated 
#'   transcript to plot.  All transcripts (assembled and annotated) overlapping
#'   \code{annotated[[ind]]} will be plotted. Default 1.
#' @param main optional character string giving the title for the resulting
#'   plot.  Default: "Assembled and Annotated Transcripts"
#' @param customCol optional vector of custom colors for the annotated
#'   transcripts. If not the same length as the number of annotated transcripts
#'   in the plot, recycling or truncation might occur.
#' 
#' @return Plots annotated transcripts on the bottom panel (shaded in gray) and
#'   assembled transcripts on the top panel (shaded with diagonal lines).
#' 
#' @author Alyssa Frazee
#' 
#' @export
#' @examples \donttest{
#' gtfPath = system.file('extdata', 'annot.gtf.gz', package='ballgown')
#' annot = gffReadGR(gtfPath, splitByTranscript=TRUE)
#' data(bg)
#' checkAssembledTx(annotated=annot, assembled=structure(bg)$trans, ind=4)
#' }
checkAssembledTx = function(assembled, annotated, ind=1, 
    main='Assembled and Annotated Transcripts', customCol=NULL){

    ol = findOverlaps(annotated, assembled)
    ol.self = findOverlaps(annotated, annotated) 
    #^ plot any overlapping transcripts also.
    indsas = subjectHits(ol)[which(queryHits(ol)==ind)]

    if(length(indsas)==0){
        warning(.makepretty('This annotated transcript was not covered by any
            assembled transcripts. No plot will be produced'), call.=FALSE)
        return(0)
    }

    # plot setup:
    par(mar=c(5,3,4,2))
    annot.list = lapply(subjectHits(ol.self)[which(queryHits(ol.self)==ind)], 
        function(i) annotated[[i]])
    names(annot.list) = names(annotated)[subjectHits(ol.self)[which(
        queryHits(ol.self)==ind)]]
    annot.df = as.data.frame(GRangesList(annot.list))
    asmbl.list = lapply(indsas, function(i) assembled[[i]])
    names(asmbl.list) = names(assembled)[indsas]
    asmbl.df = as.data.frame(GRangesList(asmbl.list))
    xax = seq(min(min(asmbl.df$start), min(annot.df$start)), 
        max(max(asmbl.df$end), max(annot.df$end)), by=1)
    if(names(asmbl.df)[2] != 'group_name'){
        stop('IRanges::as.data.frame has changed, please report')
    }
    numtx = length(unique(annot.df$group_name))+
        length(unique(asmbl.df$group_name))

    # plot base:
    plot(xax, rep(0, length(xax)), ylim=c(0.5, numtx+1), type="n", 
        xlab="genomic position", yaxt="n", ylab="")

    # plot the annotated transcripts:
    txind = 0
    for(tx in unique(annot.df$group_name)){
        txind = txind+1 #counter
        gtsub = annot.df[annot.df$group_name==tx,]
        gtsub = gtsub[order(gtsub$start),]
        for(exind in 1:dim(gtsub)[1]){
            polygon(x=c(gtsub$start[exind], gtsub$start[exind], 
                gtsub$end[exind], gtsub$end[exind]), 
                y=c(txind-0.4, txind+0.4, txind+0.4, txind-0.4), col="gray30")
            if(exind!=dim(gtsub)[1]){
                lines(c(gtsub$end[exind],gtsub$start[exind+1]), c(txind, txind),
                    lty=2, col="gray50")
            }
        }
        st.compare = gtsub$start[-1]
        en.compare = gtsub$end[-length(gtsub$end)]
        diffs = st.compare-en.compare
        if(any(diffs<0)){
            warning(paste("overlapping exons in annotated transcript",tx), 
                call.=FALSE)
        }
    }
    colorind = 0
    txind = txind+0.5
    abline(h=txind+0.25)

    # plot the assembled transcripts
    for(tx in unique(asmbl.df$group_name)){
        txind = txind+1
        gtsub = asmbl.df[asmbl.df$group_name==tx,]
        gtsub = gtsub[order(gtsub$start),]
        if(is.null(customCol)){
            mycolor = 'black'
            density = 15
        }else{
            mycolor = customCol[(colorind %% length(customCol))+1]
            density = NULL
        }
        mycolor = ifelse(is.null(customCol), 'black', 
            customCol[(colorind %% length(customCol))+1])
        for(exind in 1:dim(gtsub)[1]){
            polygon(x=c(gtsub$start[exind], gtsub$start[exind], 
                    gtsub$end[exind], gtsub$end[exind]), 
                y=c(txind-0.4, txind+0.4, txind+0.4, txind-0.4), 
                col=mycolor, density=density)
            if(exind!=dim(gtsub)[1]){
                lines(c(gtsub$end[exind],gtsub$start[exind+1]), c(txind, txind),
                    lty=2, col="gray50")
            }
        }
        st.compare = gtsub$start[-1]
        en.compare = gtsub$end[-length(gtsub$end)]
        diffs = st.compare-en.compare
        if(any(diffs<0)){
            warning(paste("overlapping exons in assembled transcript",tx), 
                call.=FALSE)
        }
        colorind = colorind+1
    }
    title(main)
    axis(side=2, at=c(median(1:length(unique(annot.df$group_name))), 
        median(1:length(unique(asmbl.df$group_name)))+
            length(unique(annot.df$group_name))+0.5), 
        labels=c("annotated", "assembled"), tick=FALSE)
}
