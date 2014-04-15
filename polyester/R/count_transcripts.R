#'<brief desc>
#'
#'<full description>
#' @param fasta <what param does>
#' @export
#' @keywords
#' @seealso
#' @return
#' @alias
#' @examples \dontrun{
#'
#'}
count_transcripts = function(fasta){
    length(readDNAStringSet(fasta))
}
