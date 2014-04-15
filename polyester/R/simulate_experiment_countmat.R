#' Simulate RNA-seq experiment 
#'
#'<full description>
#' @param fasta path to FASTA file containing transcripts from which to simulate reads 
#' @param readmat matrix with rows representing transcripts and columns representing samples. 
#' Entry i,j specifies how many reads to simulate from transcript i for sample j.
#' @param outdir character, path to folder where simulated reads should be written. Should end with "/" if specified. If unspecified, reads are written to the working directory.
#' @param fraglen Mean RNA fragment length. Sequences will be read off the end(s) of these fragments.
#' @param fragsd Standard deviation of fragment lengths. 
#' @param readlen Read length
#' @param error_rate=0.005 Sequencing error rate. Must be between 0 and 1. A uniform error model is assumed. 
#' @param paired If \code{TRUE}, paired-end reads are simulated; else single-end reads are simulated.
#' @export
#' @examples \dontrun{
#' fastapath = system.file("data", "chr22.fa", package="polyester")
#' numtx = count_transcripts(fastapath)
#' readmat = matrix(20, ncol=10, nrow=numtx)
#' readmat[1:30, 1:5] = 40
#' 
#' simulate_experiment_countmat(fastapath, readmat, outdir="./data/")
#'}
simulate_experiment_countmat = function(fasta, readmat, outdir="", 
    fraglen=250, fragsd=25, readlen=100, error_rate=0.005, paired=TRUE){

    transcripts = readDNAStringSet(fasta)

    for(i in 1:ncol(readmat)){
        tObj = rep(transcripts,times=readmat[,i])
  
        #get fragments
        tFrags = generate_fragments(tObj, fraglen=fraglen)

        #reverse_complement some of those fragments
        rctFrags = reverse_complement(tFrags)

        #add sequencing error
        errFrags = add_error(rctFrags, error_rate)

        #get reads from fragments
        reads = get_reads(errFrags, readlen, paired)

        #write read pairs
        write_reads(reads, readlen=readlen, 
            fname=paste0(outdir, 'sample_', sprintf('%02d', i)), 
            paired=paired)
    }
}