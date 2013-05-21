checkAssembledTx = function(tx.assembled, tx.annotated, ind = 1){
  #tx.assembled & tx.used should be GRangesList objects containing the assembled transcripts (this is in the ballgown object) and the annotated transcripts (user needs to make this - maybe I will write a function later).  ind tells which annotated transcript you would like to check.
  
  ol = findOverlaps(tx.annotated, tx.assembled)
  ol.self = findOverlaps(tx.annotated, tx.annotated) #plot any overlapping transcripts also.
  inds.asmb = subjectHits(ol)[which(queryHits(ol)==ind)]
  
  if(length(inds.asmb)==0){
    warning("This annotated transcript was not covered by any assembled transcripts. No plot will be produced", call.=FALSE)
    return(0)
  }
  
  # plot setup:
  par(mar=c(5,3,4,2))
  annot.list = lapply(subjectHits(ol.self)[which(queryHits(ol.self)==ind)], function(i) tx.annotated[[i]])
  names(annot.list) = names(tx.annotated)[subjectHits(ol.self)[which(queryHits(ol.self)==ind)]]
  annot.df = as.data.frame(GRangesList(annot.list))
  asmbl.list = lapply(inds.asmb, function(i) tx.assembled[[i]])
  names(asmbl.list) = names(tx.assembled)[inds.asmb]
  asmbl.df = as.data.frame(GRangesList(asmbl.list))
  xax = seq(min(min(asmbl.df$start), min(annot.df$start)), max(max(asmbl.df$end), max(annot.df$end)), by=1)
  numtx = length(unique(annot.df$element))+length(unique(asmbl.df$element))
  
  # plot base:
  plot(xax, rep(0, length(xax)), ylim=c(0.5, numtx+1), type="n", xlab="genomic position", yaxt="n", ylab="")
  
  # plot the annotated transcripts:
  txind = 0
  for(tx in unique(annot.df$element)){
    txind = txind+1 #counter
    gtsub = annot.df[annot.df$element==tx,]
    gtsub = gtsub[order(gtsub$start),]
    for(exind in 1:dim(gtsub)[1]){
      polygon(x=c(gtsub$start[exind], gtsub$start[exind], gtsub$end[exind], gtsub$end[exind]), y=c(txind-0.4, txind+0.4, txind+0.4, txind-0.4), col="gray30")
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
      polygon(x=c(gtsub$start[exind], gtsub$start[exind], gtsub$end[exind], gtsub$end[exind]), y=c(txind-0.4, txind+0.4, txind+0.4, txind-0.4), col="black", density=15, angle=45)
      if(exind!=dim(gtsub)[1]){
        lines(c(gtsub$end[exind],gtsub$start[exind+1]),c(txind, txind), lty=2, col="gray50")
      }
    }
    st.compare = gtsub$start[-1]
    en.compare = gtsub$end[-length(gtsub$end)]
    diffs = st.compare-en.compare
    if(any(diffs<0)) warning(paste("overlapping exons in assembled transcript",tx), call.=FALSE)
  }
  title("Assembled and Annotated Transcripts")
  axis(side=2, at=c(median(1:length(unique(annot.df$element))), median(1:length(unique(asmbl.df$element)))+length(unique(annot.df$element))+0.5), labels=c("annotated", "assembled"), tick=FALSE)
}
