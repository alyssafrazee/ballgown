# count the transcripts in a FASTA file
#'
#' determine how many transcripts are annotated in a FASTA file
#' @param fasta character, path to a file in FASTA format
#' @export
#' @return Number of transcripts annotated in \code{fasta}
#' @examples 
#' fastapath = system.file("data", "chr22.fa", package="polyester")
#' count_transcripts(fastapath) #918
#'
count_transcripts = function(fasta){
    length(readDNAStringSet(fasta))
}
