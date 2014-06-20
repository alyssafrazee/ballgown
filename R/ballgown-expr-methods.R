#' extract exon-level expression measurements from ballgown objects
#' 
#' @name eexpr
#' @exportMethod eexpr
#' @docType methods
#' @rdname eexpr
#' @aliases eexpr,ballgown-method
#' @param x a ballgown object
#' @param meas type of measurement to extract. Can be "rcount", "ucount", "mrcount", "cov", "mcov",
#'             or "all". Default "rcount".
#' @return exon-by-sample matrix containing exon-level expression values (measured by \code{meas})
setMethod('eexpr', 'ballgown', function(x, meas='rcount'){
    emeas = c('all', 'rcount', 'ucount', 'mrcount', 'cov', 'mcov', 'cov_sd', 'mcov_sd')
    meas = match.arg(meas, emeas)
    if(meas!='all'){
        if(!identical(x@meas, 'all')){
            if(!(meas %in% x@meas)){
                meas_avail = paste(c(intersect(x@meas, emeas), 'all'), collapse=', ')
                msg = paste(meas, 'measurements were not included when this ballgown object was
                    created. Instead, you can use eexpr with these measurements:', meas_avail)
                stop(.makepretty(msg))
            }
        }
        expr = data(x)$exon[,-c(1:5)]
        expr = expr[,sapply(colnames(expr), function(x) strsplit(x,split="\\.")[[1]][1]==meas)]
        rownames(expr) = data(x)$exon$e_id
        expr = as.matrix(expr)
    }else{
        if(!any(emeas[-1] %in% x@meas) & !identical(x@meas, 'all')){
            msg = 'No exon-level expression measurements are included in this ballgown object. 
                Returning data frame of all exons\' genomic coordinates.'
            warning(.makepretty(msg))
        }
        expr = data(x)$exon
    }
    return(expr)
})

#' extract transcript-level expression measurements from ballgown objects
#' 
#' @name texpr
#' @exportMethod texpr
#' @docType methods
#' @rdname texpr
#' @aliases texpr,ballgown-method
#' @param x a ballgown object
#' @param meas type of measurement to extract. Can be "cov", "FPKM", or "all". Default "FPKM".
#' @return transcript-by-sample matrix containing expression values (measured by \code{meas})
setMethod('texpr', 'ballgown', function(x, meas='FPKM'){
    meas = match.arg(meas, c('cov', 'FPKM', 'all'))
    if(meas!='all'){
        if(!identical(x@meas, 'all')){
            if(!(meas %in% x@meas)){
                meas_avail = paste(c(intersect(x@meas, c('cov', 'FPKM')), 'all'), collapse=', ')
                msg = paste(meas, 'measurements were not included when this ballgown object was
                    created. Instead, you can use texpr with these measurements:', meas_avail)
                stop(.makepretty(msg))
            }
        }
        expr = data(x)$trans[,-c(1:10)]
        expr = expr[,sapply(colnames(expr), function(x) strsplit(x,split="\\.")[[1]][1]==meas)]
        rownames(expr) = data(x)$trans$t_id
        expr = as.matrix(expr)
    }else{
        if(!any(c('cov', 'FPKM') %in% x@meas) & !identical(x@meas, 'all')){
            msg = 'No transcript-level expression measurements are included in this ballgown 
                object. Returning data frame of transcript structure information.'
            warning(.makepretty(msg))
        }
        expr = data(x)$trans
    }
    return(expr)
})

#' extract transcript-level expression measurements from ballgown objects
#' 
#' @name iexpr
#' @exportMethod iexpr
#' @docType methods
#' @rdname iexpr
#' @aliases iexpr,ballgown-method
#' @param x a ballgown object
#' @param meas type of measurement to extract. Can be "rcount", "ucount", "mrcount", or "all". 
#'   Default "rcount".
#' @return intron-by-sample matrix containing the number of reads (measured as specified by 
#'   \code{meas}) supporting each intron, in each sample.
setMethod('iexpr', 'ballgown', function(x, meas='rcount'){
    meas = match.arg(meas, c('rcount', 'ucount', 'mrcount', 'all'))
    if(meas!='all'){
        if(!identical(x@meas, 'all')){
            if(!(meas %in% x@meas)){
                meas_avail = paste(c(intersect(x@meas, c('rcount', 'ucount', 'mrcount')), 'all'),
                    collapse=', ')
                msg = paste(meas, 'measurements were not included when this ballgown object was
                    created. Instead, you can use iexpr with these measurements:', meas_avail)
                stop(.makepretty(msg))
            }
        }
        expr = data(x)$intron[,-c(1:5)]
        expr = expr[,sapply(colnames(expr), function(x) strsplit(x,split="\\.")[[1]][1]==meas)]
        rownames(expr) = data(x)$intron$i_id
        expr = as.matrix(expr)
    }else{
        if(!any(c('ucount', 'rcount', 'mrcount') %in% x@meas) & !identical(x@meas, 'all')){
            msg = 'No intron-level expression measurements are included in this ballgown object. 
                Returning data frame of all introns\' genomic coordinates.'
            warning(.makepretty(msg))
        }
        expr = data(x)$intron
    }
    return(expr)
})

#' extract gene-level expression measurements from ballgown objects
#' 
#' gene-level measurements are done by appropriately combining FPKMs from the transcripts comprising
#' the gene. 
#'
#' @name gexpr
#' @exportMethod gexpr
#' @docType methods
#' @rdname gexpr
#' @aliases gexpr,ballgown-method
#' @param x a ballgown object
#' @return gene-by-sample matrix containing per-sample gene FPKMs.
setMethod('gexpr', 'ballgown', function(x){
    if(!('FPKM' %in% x@meas) & !identical(x@meas, 'all')){
        msg = 'gene expression measurements can only be calculated for ballgown objects including 
            transcript-level FPKM values.' 
        stop(.makepretty(msg))
    }
    gnames = indexes(x)$t2g$g_id
    inds_by_gene = split(seq(along=gnames), gnames)
    tmeas = texpr(x, "FPKM")
    gid_by_exon = lapply(1:nrow(tmeas), function(i){
        rep(texpr(x, 'all')$gene_id[i], texpr(x, 'all')$num_exons[i])})
    ulstruct = unlist(structure(x)$trans)
    glist = split(ulstruct, unlist(gid_by_exon))
    glengths = sapply(width(reduce(glist)), sum)
    tlengths = sapply(width(structure(x)$trans), sum)
    tfrags = tlengths/1000 * tmeas
    expr = t(sapply(1:length(inds_by_gene), function(i){
        colSums(tfrags[inds_by_gene[[i]],,drop=FALSE]) / glengths[i]}))
    rownames(expr) = names(inds_by_gene)
    return(expr)
})


