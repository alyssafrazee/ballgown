#'<brief desc>
#'
#'<full description>
#' @param fasta <what param does>
#' @export
#' @return Number of transcripts annotated in \code{fasta}
#' @examples \dontrun{
#'
#'}
count_transcripts = function(fasta){
    length(readDNAStringSet(fasta))
}
