#' write sequencing reads to disk
#'
#' given a DNAStringSet representing simulated sequencing reads, write FASTA files to disk
#' representing the simulated reads.
#' @param reads DNAStringSet representing sequencing reads
#' @param fname file path to where sequencing reads should be written. Should be not contain 
#' ".fasta" (this is appended automatically). 
#' @param readlen maximum length of the reads in \code{reads}.
#' @param paired If \code{TRUE}, reads are assumed to be in pairs: i.e., read 1 and read 2 in 
#' \code{reads} are the left and right mate (respectively) of a read pair; same with read 3 and 
#' read 4, etc. The odd-numbered reads are written to \code{fname_1.fasta} and the even-numbered
#' reads are written to \code{fname_2.fasta}. If \code{FALSE}, reads are assumed to be single-end
#' and just one file, \code{fname.fasta}, is written.
#' @export
#' @details The \code{\link{get_reads}} function returns a DNAStringSet object representing 
#' sequencing reads that can be directly passed to \code{write_reads}. If output other than that 
#' from \code{get_reads} is used and \code{paired} is \code{TRUE}, make sure \code{reads} is 
#' ordered properly (i.e., that mate pairs appear together and that the left mate appears first).
#' @seealso \code{\link{get_reads}}
write_reads = function(reads, fname, readlen, paired = TRUE){
    if(paired){
        lefts = seq(1, length(reads), by=2)
        rights = seq(2, length(reads), by=2)
        outfl = paste0(fname, '_1.fasta')
        outfr = paste0(fname, '_2.fasta')
        writeXStringSet(reads[lefts], file=outfl, format="fasta", width=readlen)
        writeXStringSet(reads[rights], file=outfr, format="fasta", width=readlen)
    }else{
        outf = paste0(fname, '.fasta')
        writeXStringSet(reads, file=outf, format="fasta", width=readlen)
    }
}
