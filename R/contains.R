#' determine if one set of GRanges fully contains any of another set of GRanges
#' 
#' @param transcripts \code{GRangesList} object (assume for now that it
#'    represents transcripts)
#' @param cds \code{GRangesList} object (assume for now that it represents sets
#'    of coding sequences)
#' 
#' @details If \code{gown} is a \code{ballgown} object, \code{transcripts} can
#'    be \code{structure(gown)$trans} (or any subset). 
#' 
#' @return vector with length equal to \code{length(transcripts)}, where each
#'    entry is \code{TRUE} if the corresponding transcript contains a coding
#'    sequence (i.e., is a superset of at least one entry of \code{cds}). 
#' 
#' @author Alyssa Frazee
#' 
#' @export
#' 
#' @examples
#' ## pretend this annotation is coding sequence:
#' gtfPath = system.file('extdata', 'annot.gtf.gz', package='ballgown')
#' annot = gffReadGR(gtfPath, splitByTranscript=TRUE)
#' data(bg)
#' results = contains(structure(bg)$trans, annot)
#' # results is a boolean vector
#' sum(results) #61
contains = function(transcripts, cds){
    
    stopifnot(class(transcripts) == 'GRangesList' & class(cds) == 'GRangesList')

    # things will get out of order: keep track
    names_transcripts = names(transcripts)
    names_cds = names(cds)

    # get rid of metadata:
    transcriptsUL = unlist(transcripts)
    mcols(transcriptsUL) = NULL
    transcripts = split(transcriptsUL, names(transcriptsUL))
    transcripts = transcripts[match(names_transcripts, names(transcripts))]

    # for assembly: get rid of metadata
    cdsUL = unlist(cds)
    mcols(cdsUL) = NULL
    cds = split(cdsUL, names(cdsUL))
    cds = cds[match(names_cds, names(cds))]

    # find the overlaps:
    ol = findOverlaps(transcripts, cds)

    # get lists of corresponding assembled/annotated transcripts:
    transcripts_sort = transcripts[queryHits(ol)]
    cds_sort = cds[subjectHits(ol)]

    # concatenate the appropriate GRanges objects (for overlapping transcripts)
    transcripts_overlapIDs = rep(1:length(ol), 
        times=elementNROWS(transcripts_sort))
    cds_overlapIDs = rep(1:length(ol), times=elementNROWS(cds_sort))
    all_transcripts = c(unlist(transcripts_sort), unlist(cds_sort))
    overlapping = split(all_transcripts, 
        c(transcripts_overlapIDs, cds_overlapIDs))
    
    # decide on containment:
    stopifnot(all(as.numeric(names(overlapping)) == 1:length(ol)))
    # names(overlapping) = 1:length(ol), 
    # corresponds to overlap IDs (confirmed in above check) 
    coverages = coverage(ranges(overlapping))
    runvals = runValue(coverages)
    runlengths = runLength(coverages)
    cds_lengths = sapply(width(cds_sort), sum)
    containsCDS = sum(runlengths[runvals==2]) == cds_lengths

    # return data in right order:
    allOL = split(containsCDS, queryHits(ol))
    ret = rep(FALSE, length(transcripts))
    ret[as.numeric(names(allOL))] = sapply(allOL, any)
    return(ret)
}

