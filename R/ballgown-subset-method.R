# subset method for ballgown
setMethod("subset", "ballgown", function(x, cond, global=TRUE){

    ctext = ifelse(global, deparse(substitute(cond)), cond) # this means that inside another function, you can make a string argument to give to subset.
    trans = subset(data(x)$trans, eval(parse(text=ctext)))  
    
    thetx = trans$t_id
    
    inttmp = split(indexes(x)$i2t$i_id, indexes(x)$i2t$t_id)
    theint = as.numeric(unique(unlist(inttmp[names(inttmp) %in% thetx])))
    intron = subset(data(x)$intron, i_id %in% theint)
    
    extmp = split(indexes(x)$e2t$e_id, indexes(x)$e2t$t_id)
    theex = as.numeric(unique(unlist(extmp[names(extmp) %in% thetx])))
    exon = subset(data(x)$exon, e_id %in% theex)
    
    e2t = subset(indexes(x)$e2t, t_id %in% thetx)
    i2t = subset(indexes(x)$i2t, t_id %in% thetx)
    t2g = subset(indexes(x)$t2g, t_id %in% thetx)
    
    introngr = structure(x)$intron[elementMetadata(structure(x)$intron)$id %in% theint]
    exongr = structure(x)$exon[elementMetadata(structure(x)$exon)$id %in% theex]
    grltxids = sapply(names(structure(x)$trans), function(a) as.numeric(substr(a, 3, nchar(a))))
    transgrl = structure(x)$trans[grltxids %in% thetx]
    
    return(new("ballgown", data = list(intron=intron, exon=exon, trans=trans), indexes = list(e2t=e2t, i2t=i2t, t2g=t2g, bamfiles = indexes(x)$bamfiles, pData = indexes(x)$pData), structure=list(intron = introngr, exon=exongr, trans=transgrl), dirs = dirs(x), mergedDate=mergedDate(x)))
} )
