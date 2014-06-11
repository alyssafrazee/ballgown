# subset method for ballgown
setMethod("subset", "ballgown", function(x, cond, genomesubset=TRUE){
    stopifnot(class(cond) == 'character')

    # if you are subsetting by something in the genome (say, a chromosome):
    if(genomesubset){
        trans = subset(data(x)$trans, eval(parse(text=cond)))  
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
        grltxids = as.numeric(names(structure(x)$trans))
        transgrl = structure(x)$trans[grltxids %in% thetx]
    
        return(new("ballgown", data=list(intron=intron, exon=exon, trans=trans), 
            indexes=list(e2t=e2t, i2t=i2t, t2g=t2g, bamfiles=indexes(x)$bamfiles, 
                pData=indexes(x)$pData), 
                structure=list(intron=introngr, exon=exongr, trans=transgrl), 
                dirs=dirs(x), mergedDate=mergedDate(x)))
    }else{
        # you're doing a phenotype subset
        # structure, some indexes, dirs, and mergedDate stay the same
        # change: data, indexes(pData), and indexes(bamfiles)
        
        ## pData
        newpd = subset(pData(x), eval(parse(text=cond)))
        
        ## bamfiles
        newsampnames = newpd[,1]
        rowIndsToKeep = which(pData(x)[,1] %in% newsampnames)
        if(!is.null(indexes(x)$bamfiles)){
            newbamfiles = indexes(x)$bamfiles[rowIndsToKeep]
        }else{
            newbamfiles = NULL
        }

        ## transcript data
        txcolsamples = sapply(names(texpr(x, 'all')), getsamp, USE.NAMES=FALSE)
        txKeepCols = c(1:10, which(txcolsamples %in% newsampnames))
        newtdat = texpr(x, 'all')[,txKeepCols]

        ## exon data
        excolsamples = sapply(names(eexpr(x, 'all')), getsamp, USE.NAMES=FALSE)
        exKeepCols = c(1:5, which(excolsamples %in% newsampnames))
        newedat = eexpr(x, 'all')[,exKeepCols]

        ## intron data
        icolsamples = sapply(names(iexpr(x, 'all')), getsamp, USE.NAMES=FALSE)
        iKeepCols = c(1:5, which(icolsamples %in% newsampnames))
        newidat = iexpr(x, 'all')[,iKeepCols]

        return(new("ballgown", data=list(intron=newidat, exon=newedat, trans=newtdat), 
            indexes=list(e2t=indexes(x)$e2t, i2t=indexes(x)$i2t, t2g=indexes(x)$t2g, 
            bamfiles=newbamfiles, pData=newpd), 
            structure=list(intron=structure(x)$intron, exon=structure(x)$exon, 
                trans=structure(x)$trans),
            dirs=dirs(x)[rowIndsToKeep], mergedDate=mergedDate(x)))
    }
} )
