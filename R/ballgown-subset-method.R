## function that returns 2nd half of expression column names
## (regardless of whether they have '.' in them)
half2 = function(x) paste(x[-1], collapse='.')


#' subset ballgown objects to specific samples or genomic locations
#' 
#' @name subset
#' @exportMethod subset
#' @docType methods
#' @rdname subset
#' @aliases subset,ballgown-method
#' @param x a ballgown object
#' @param cond Condition on which to subset. See details.
#' @param genomesubset if TRUE, subset \code{x} to a specific part of the 
#'   genome. Otherwise, subset x to only include specific samples. TRUE by 
#'   default.
#' @return a subsetted ballgown object, containing only the regions or samples
#'   satisfying \code{cond}.
#' 
#' @details To use \code{subset}, you must provide the \code{cond} argument as a
#'   string representing a logical expression specifying your desired subset.
#'   The subset expression can either involve column names of
#'   \code{texpr(x, "all")} (if \code{genomesubset} is \code{TRUE}) or of 
#'   \code{pData(x)} (if \code{genomesubset} is \code{FALSE}). For example, if
#'   you wanted a ballgown object for only chromosome 22, you might call 
#'   \code{subset(x, "chr == 'chr22'")}. (Be sure to handle quotes within
#'   character strings appropriately). 
#' 
#' @author Alyssa Frazee
#' 
#' @examples
#' data(bg)
#' bg_twogenes = subset(bg, "gene_id=='XLOC_000454' | gene_id=='XLOC_000024'")
#' bg_twogenes 
#' # ballgown instance with 4 assembled transcripts and 20 samples
#' 
#' bg_group0 = subset(bg, "group == 0", genomesubset=FALSE)
#' bg_group0 
#' # ballgown instance with 100 assembled transcripts and 10 samples
setMethod("subset", "ballgown", function(x, cond, genomesubset=TRUE){
    stopifnot(class(cond) == 'character')

    # if you are subsetting by something in the genome (say, a chromosome):
    if(genomesubset){
        trans = subset(expr(x)$trans, eval(parse(text=cond)))  
        thetx = trans$t_id
    
        inttmp = split(indexes(x)$i2t$i_id, indexes(x)$i2t$t_id)
        theint = as.numeric(unique(unlist(inttmp[names(inttmp) %in% thetx])))
        intron = subset(expr(x)$intron, i_id %in% theint)
    
        extmp = split(indexes(x)$e2t$e_id, indexes(x)$e2t$t_id)
        theex = as.numeric(unique(unlist(extmp[names(extmp) %in% thetx])))
        exon = subset(expr(x)$exon, e_id %in% theex)
    
        e2t = subset(indexes(x)$e2t, t_id %in% thetx)
        i2t = subset(indexes(x)$i2t, t_id %in% thetx)
        t2g = subset(indexes(x)$t2g, t_id %in% thetx)
    
        introngr = structure(x)$intron[elementMetadata(structure(x)$intron)$id 
            %in% theint]
        exongr = structure(x)$exon[elementMetadata(structure(x)$exon)$id 
            %in% theex]
        grltxids = as.numeric(names(structure(x)$trans))
        transgrl = structure(x)$trans[grltxids %in% thetx]
    
        return(new("ballgown", expr=list(intron=intron, exon=exon, trans=trans),
            indexes=list(e2t=e2t, i2t=i2t, t2g=t2g, 
                bamfiles=indexes(x)$bamfiles, pData=indexes(x)$pData), 
            structure=list(intron=introngr, exon=exongr, trans=transgrl), 
            dirs=dirs(x), mergedDate=mergedDate(x), meas=x@meas, RSEM=x@RSEM))
    }else{
        # you're doing a phenotype subset
        # structure, some indexes, dirs, and mergedDate stay the same
        # change: data, indexes(pData), and indexes(bamfiles)
        
        ## pData
        newpd = subset(pData(x), eval(parse(text=cond)))
        newpd = droplevels(newpd)
        
        ## bamfiles
        newsampnames = newpd[,1]
        rowIndsToKeep = which(pData(x)[,1] %in% newsampnames)
        if(!is.null(indexes(x)$bamfiles)){
            newbamfiles = indexes(x)$bamfiles[rowIndsToKeep]
        }else{
            newbamfiles = NULL
        }

        ## transcript data
        txcolsamples = sapply(strsplit(names(texpr(x, 'all')), '\\.'), half2)
        txKeepCols = c(1:10, which(txcolsamples %in% newsampnames))
        newtdat = texpr(x, 'all')[,txKeepCols]

        ## exon data
        excolsamples = sapply(strsplit(names(eexpr(x, 'all')), '\\.'), half2)
        exKeepCols = c(1:5, which(excolsamples %in% newsampnames))
        newedat = eexpr(x, 'all')[,exKeepCols]

        ## intron data
        icolsamples = sapply(strsplit(names(iexpr(x, 'all')), '\\.'), half2)
        iKeepCols = c(1:5, which(icolsamples %in% newsampnames))
        newidat = iexpr(x, 'all')[,iKeepCols]

        return(new("ballgown", 
            expr=list(intron=newidat, exon=newedat, trans=newtdat), 
            indexes=list(e2t=indexes(x)$e2t, i2t=indexes(x)$i2t,
                t2g=indexes(x)$t2g, bamfiles=newbamfiles, pData=newpd), 
            structure=list(intron=structure(x)$intron, exon=structure(x)$exon, 
                trans=structure(x)$trans),
            dirs=dirs(x)[rowIndsToKeep], mergedDate=mergedDate(x), 
            meas=x@meas, RSEM=x@RSEM))
    }
} )
