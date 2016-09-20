#' match assembled transcripts to annotated transcripts
#'
#' @param assembled \code{GRangesList} object representing assembled transcripts
#' @param annotated \code{GRangesList} object representing annotated transcripts
#' 
#' @details If \code{gown} is a \code{ballgown} object, \code{assembled} can be 
#'   \code{structure(gown)$trans} (or any subset). You can generate a 
#'   \code{GRangesList} object 
#' containing annotated transcripts from a gtf file using the 
#'   \code{\link{gffReadGR}} function and 
#' setting \code{splitByTranscripts=TRUE}.
#' 
#' @return data frame, where each row contains \code{assembledInd} and 
#'   \code{annotatedInd} (indexes of overlapping transcripts in \code{assembled}
#'   and \code{annotated}), and the percent overlap between the two transcripts.
#' 
#' @author Alyssa Frazee
#' 
#' @export
#' @examples
#' data(bg)
#' gtfPath = system.file('extdata', 'annot.gtf.gz', package='ballgown')
#' annot = gffReadGR(gtfPath, splitByTranscript=TRUE)
#' info = annotate_assembly(assembled=structure(bg)$trans, annotated=annot)

annotate_assembly = function(assembled, annotated){
    
    stopifnot(class(annotated) == 'GRangesList' & 
        class(assembled) == 'GRangesList')

    # things will get out of order: keep track
    names_annot = names(annotated)
    names_assmb = names(assembled)

    # get rid of metadata:
    annotatedUL = unlist(annotated)
    mcols(annotatedUL) = NULL
    annotated = split(annotatedUL, names(annotatedUL))
    annotated = annotated[match(names_annot, names(annotated))]

    # for assembly: get rid of metadata
    assembledUL = unlist(assembled)
    mcols(assembledUL) = NULL
    assembled = split(assembledUL, names(assembledUL))
    assembled = assembled[match(names_assmb, names(assembled))]

    # find the overlaps:
    ol = findOverlaps(assembled, annotated)
    ol = ol[order(queryHits(ol), subjectHits(ol))]

    # get lists of corresponding assembled/annotated transcripts:
    assembled_sort = assembled[queryHits(ol)]
    annotated_sort = annotated[subjectHits(ol)]

    # concatenate the appropriate GRanges objects (for overlapping transcripts)
    assembled_overlapIDs = rep(1:length(ol), 
        times=elementNROWS(assembled_sort))
    annotated_overlapIDs = rep(1:length(ol), 
        times=elementNROWS(annotated_sort))
    all_transcripts = c(unlist(assembled_sort), unlist(annotated_sort))
    overlapping = split(all_transcripts, 
        c(assembled_overlapIDs, annotated_overlapIDs))
    
    # calculate percentage overlap:
    coverages = coverage(ranges(overlapping))
    runvals = runValue(coverages)
    runlengths = runLength(coverages)
    pct = sum(runlengths[runvals==2]) / sum(runlengths[runvals==1 | runvals==2])

    return(data.frame(assembledInd=queryHits(ol), 
        annotatedInd=subjectHits(ol), percent=pct))
}

