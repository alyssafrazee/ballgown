#' @title subset ballgown objects using an expression filter
#' 
#' @description Create a new ballgown object containing only transcripts passing
#'   a mean expression filter
#' @param gown a ballgown object
#' @param cutoff transcripts must have mean expression across samples above this
#'   value to be included in the return
#' @param meas how should transcript expression be measured? Default FPKM, but 
#'   can also be \code{'cov'}.
#' @export
#' @seealso \code{\link{subset}}
#' @return A new ballgown object derived from \code{gown}, but only containing 
#'   transcripts (and associated exons/introns) with mean \code{meas} greater
#'   than \code{cutoff} across all samples.
#' @examples
#'   data(bg)
#'   # make a ballgown object containing only transcripts with mean FPKM > 100:
#'   over100 = exprfilter(bg, cutoff=100)  
#'
exprfilter = function(gown, cutoff, meas='FPKM'){
    e_id = t_id = i_id = NULL #prep for use in subset function / R CMD CHECK
    meas = match.arg(meas, c('FPKM', 'cov'))

    stopifnot(gown@meas == 'all' | meas %in% gown@meas)

    means = rowMeans(texpr(gown, meas))
    hiexpr_index = which(means > cutoff)
    trans = texpr(gown, 'all')[hiexpr_index,]
    thetx = trans$t_id

    inttmp = split(indexes(gown)$i2t$i_id, indexes(gown)$i2t$t_id)
    theint = as.numeric(unique(unlist(inttmp[names(inttmp) %in% thetx])))
    intron = subset(expr(gown)$intron, i_id %in% theint)

    extmp = split(indexes(gown)$e2t$e_id, indexes(gown)$e2t$t_id)
    theex = as.numeric(unique(unlist(extmp[names(extmp) %in% thetx])))
    exon = subset(expr(gown)$exon, e_id %in% theex)

    e2t = subset(indexes(gown)$e2t, t_id %in% thetx)
    i2t = subset(indexes(gown)$i2t, t_id %in% thetx)
    t2g = subset(indexes(gown)$t2g, t_id %in% thetx)

    introngr = structure(gown)$intron[elementMetadata(structure(gown)$intron)$id
        %in% theint]
    exongr = structure(gown)$exon[elementMetadata(structure(gown)$exon)$id 
        %in% theex]
    transgrl = structure(gown)$trans[names(structure(gown)$trans) %in% thetx]
    return(new('ballgown', 
        expr=list(intron=intron, exon=exon, trans=trans), 
        indexes=list(e2t=e2t, i2t=i2t, t2g=t2g, 
            bamfiles=indexes(gown)$bamfiles, pData=indexes(gown)$pData), 
        structure=list(intron=introngr, exon=exongr, trans=transgrl), 
        dirs=dirs(gown), mergedDate=mergedDate(gown), RSEM=gown@RSEM, 
        meas=gown@meas))
}

