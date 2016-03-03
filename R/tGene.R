#' @title Connect a transcript to its gene
#'
#' @description find the gene to which a transcript belongs
#' @param bg ballgown object
#' @param transcript transcript identifier
#' @param tid set to \code{TRUE} if \code{transcript} is a numeric transcript
#'   identifier (i.e., \code{t_id} in expression tables), or \code{FALSE} if 
#'   \code{transcript} is a named identifie (e.g., \code{TCONS_000001} or 
#'   similar.
#' @param gid if \code{FALSE}, return the gene *name* associated with 
#'   \code{transcript} in \code{bg} instead of the gene *id*, which is returned
#'   by default. Take care to remember that not all ballgown objects include 
#'   gene *name* information. (They do all include gene IDs). 
#' @param warnme if \code{TRUE}, and if \code{gid} is \code{FALSE}, print a 
#'   warning if no gene name is available for the transcript. This could either
#'   mean the transcript didn't overlap an annotated gene, or that no gene names
#'   were included when \code{bg} was created.
#' @export
#' @examples
#'   data(bg)
#'   tGene(bg, 10)
#'   tGene(bg, 'TCONS_00000010', tid=FALSE)
#'   tGene(bg, 10, gid=FALSE) #empty: no gene names included in bg.
#'

tGene = function(bg, transcript, tid=TRUE, gid=TRUE, warnme=TRUE){
    if(!(transcript %in% transcriptNames(bg)) &
        !(transcript %in% transcriptIDs(bg))){
        stop('transcript is not contained in ballgown object.')
    }

    expr_col_ind = ifelse(tid, which(colnames(texpr(bg,'all')) == 't_id'),
        which(colnames(texpr(bg,'all')) == 't_name'))
    i = which(texpr(bg, 'all')[,expr_col_ind] == transcript)

    if(gid){
        return(texpr(bg, 'all')$gene_id[i])
    }else{
        ret = texpr(bg, 'all')$gene_name[i]
        if(ret=="" & warnme){
            warning(.makepretty('No gene name available for this transcript.
                To return gene ID, set gid=TRUE.'))
        }
        return(ret)
    }
}
