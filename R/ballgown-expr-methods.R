#' extract exon-level expression measurements from ballgown objects
#' 
#' @name eexpr
#' @exportMethod eexpr
#' @docType methods
#' @rdname eexpr
#' @aliases eexpr,ballgown-method
#' @param x a ballgown object
#' @param meas type of measurement to extract. Can be "rcount", "ucount", 
#'   "mrcount", "cov", "mcov", or "all". Default "rcount".
#' @return exon-by-sample matrix containing exon-level expression values
#'   (measured by \code{meas}). If \code{meas} is \code{"all"}, or 
#'   \code{x@@RSEM} is TRUE, a data frame is returned, containing all
#'   measurements and location information.
#' 
#' @examples 
#' data(bg)
#' exon_rcount_matrix = eexpr(bg)
#' exon_ucount_matrix = eexpr(bg, 'ucount')
#' exon_data_frame = eexpr(bg, 'all')
setMethod('eexpr', 'ballgown', function(x, meas='rcount'){
    if(x@RSEM){
        return(expr(x)$exon)
    }
    emeas = c('all', 'rcount', 'ucount', 'mrcount', 'cov', 'mcov', 'cov_sd', 
        'mcov_sd')
    meas = match.arg(meas, emeas)
    if(meas!='all'){
        if(!identical(x@meas, 'all')){
            if(!(meas %in% x@meas)){
                meas_avail = paste(c(intersect(x@meas, emeas), 'all'), 
                    collapse=', ')
                msg = paste(meas, 'measurements were not included when this
                    ballgown object was created. Instead, you can use eexpr
                    with these measurements:', meas_avail)
                stop(.makepretty(msg))
            }
        }
        mat = expr(x)$exon[,-c(1:5)]
        mat = subset(mat, select=paste(meas, sampleNames(x), sep='.'))
        rownames(mat) = expr(x)$exon$e_id
        mat = as.matrix(mat)
    }else{
        if(!any(emeas[-1] %in% x@meas) & !identical(x@meas, 'all')){
            msg = 'No exon-level expression measurements are included in this
                    ballgown object. Returning data frame of all exons\' genomic
                    coordinates.'
            warning(.makepretty(msg))
        }
        mat = expr(x)$exon
    }
    return(mat)
})

#' extract transcript-level expression measurements from ballgown objects
#' 
#' @name texpr
#' @exportMethod texpr
#' @docType methods
#' @rdname texpr
#' @aliases texpr,ballgown-method
#' @param x a ballgown object
#' @param meas type of measurement to extract. Can be "cov", "FPKM", or "all". 
#'   Default "FPKM".
#' @return transcript-by-sample matrix containing expression values (measured by
#'   \code{meas}). If \code{meas} is \code{"all"}, a data frame is returned, 
#'   containing all measurements and location information.
#' 
#' @examples
#' data(bg)
#' transcript_fpkm_matrix = texpr(bg)
#' transcript_data_frame = texpr(bg, 'all')
setMethod('texpr', 'ballgown', function(x, meas='FPKM'){
    if(x@RSEM){
        meas = match.arg(meas, c('TPM', 'FPKM', 'all'))
    }else{
        meas = match.arg(meas, c('cov', 'FPKM', 'all'))
    }
    if(meas!='all'){
        if(!identical(x@meas, 'all')){
            if(!(meas %in% x@meas)){
                meas_avail = paste(c(intersect(x@meas, c('cov', 'FPKM', 'TPM')),
                    'all'), collapse=', ')
                msg = paste(meas, 'measurements were not included when this
                    ballgown object was created. Instead, you can use texpr with
                    these measurements:', meas_avail)
                stop(.makepretty(msg))
            }
        }
        mat = expr(x)$trans[,-c(1:10)]
        mat = subset(mat, select=paste(meas, sampleNames(x), sep='.'))
        rownames(mat) = expr(x)$trans$t_id
        mat = as.matrix(mat)
    }else{
        if(!any(c('cov', 'FPKM', 'TPM') %in% x@meas) & 
                !identical(x@meas, 'all')){
            msg = 'No transcript-level expression measurements are included in
            this ballgown object. Returning data frame of transcript structure
            information.'
            warning(.makepretty(msg))
        }
        mat = expr(x)$trans
    }
    return(mat)
})

#' extract transcript-level expression measurements from ballgown objects
#' 
#' @name iexpr
#' @exportMethod iexpr
#' @docType methods
#' @rdname iexpr
#' @aliases iexpr,ballgown-method
#' @param x a ballgown object
#' @param meas type of measurement to extract. Can be "rcount", "ucount", 
#'   "mrcount", or "all". Default "rcount".
#' @return intron-by-sample matrix containing the number of reads (measured as
#'   specified by \code{meas}) supporting each intron, in each sample. If
#'   \code{meas} is \code{"all"}, a data frame is returned, containing all
#'   measurements and location information.
#' 
#' @examples
#' data(bg)
#' intron_rcount_matrix = iexpr(bg)
#' intron_data_frame = iexpr(bg, 'all')
setMethod('iexpr', 'ballgown', function(x, meas='rcount'){
    if(x@RSEM){
        return(expr(x)$intron)
    }
    meas = match.arg(meas, c('rcount', 'ucount', 'mrcount', 'all'))
    if(meas!='all'){
        if(!identical(x@meas, 'all')){
            if(!(meas %in% x@meas)){
                meas_avail = paste(c(intersect(x@meas, 
                    c('rcount', 'ucount', 'mrcount')), 'all'), collapse=', ')
                msg = paste(meas, 'measurements were not included when this
                    ballgown object was created. Instead, you can use iexpr with
                    these measurements:', meas_avail)
                stop(.makepretty(msg))
            }
        }
        mat = expr(x)$intron[,-c(1:5)]
        mat = subset(mat, select=paste(meas, sampleNames(x), sep='.'))
        rownames(mat) = expr(x)$intron$i_id
        mat = as.matrix(mat)
    }else{
        if(!any(c('ucount', 'rcount', 'mrcount') %in% x@meas) & 
            !identical(x@meas, 'all')){
            msg = 'No intron-level expression measurements are included in this
            ballgown object. Returning data frame of all introns\' genomic
            coordinates.'
            warning(.makepretty(msg))
        }
        mat = expr(x)$intron
    }
    return(mat)
})

#' extract gene-level expression measurements from ballgown objects
#' 
#' For objects created with Cufflinks/Tablemaker, gene-level measurements are 
#'   calculated by appropriately combining FPKMs from the transcripts comprising
#'   the gene. For objects created with RSEM, gene-level measurements are
#'   extracted directly from the RSEM output.
#'
#' @name gexpr
#' @exportMethod gexpr
#' @docType methods
#' @rdname gexpr
#' @aliases gexpr,ballgown-method
#' @param x a ballgown object
#' @return gene-by-sample matrix containing per-sample gene measurements.
#' 
#' @examples
#' data(bg)
#' gene_matrix = gexpr(bg)
setMethod('gexpr', 'ballgown', function(x){
    if(x@RSEM){
        return(expr(x)$gm)
    }else{
        if(!('FPKM' %in% x@meas) & !identical(x@meas, 'all')){
            msg = 'gene expression measurements can only be calculated for
            ballgown objects including transcript-level FPKM values.' 
            stop(.makepretty(msg))
        }
        gnames = as.character(indexes(x)$t2g$g_id)
        inds_by_gene = split(seq(along=gnames), gnames)
        tmeas = texpr(x, "FPKM")
        gid_by_exon = lapply(1:nrow(tmeas), function(i){
            rep(texpr(x, 'all')$gene_id[i], texpr(x, 'all')$num_exons[i])})
        ulstruct = unlist(structure(x)$trans)
        glist = split(ulstruct, unlist(gid_by_exon))
        glengths = sapply(width(reduce(glist)), sum)
        tlengths = sapply(width(structure(x)$trans), sum)
        tfrags = tlengths * tmeas
        mat = lapply(1:length(inds_by_gene), function(i){
            colSums(tfrags[inds_by_gene[[i]],,drop=FALSE]) / glengths[i]})
        mat = do.call(rbind, mat)
        rownames(mat) = names(inds_by_gene)
        return(mat)
    }
})


