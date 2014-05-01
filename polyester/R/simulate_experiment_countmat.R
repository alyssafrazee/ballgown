#' Simulate RNA-seq experiment 
#'
#' create FASTA files containing RNA-seq reads simulated from provided transcripts, with optional 
#' differential expression between two groups (designated via read count matrix)
#' @param fasta path to FASTA file containing transcripts from which to simulate reads. See details.
#' @param gtf path to GTF file containing transcript structures from which reads should be 
#' simulated. See details.
#' @param seqpath path to folder containing one FASTA file (\code{.fa} extension) for each 
#' chromosome in \code{gtf}. See details. 
#' @param readmat matrix with rows representing transcripts and columns representing samples. 
#' Entry i,j specifies how many reads to simulate from transcript i for sample j.
#' @param outdir character, path to folder where simulated reads should be written. Should end with 
#' "/" if specified. If unspecified, reads are written to the working directory.
#' @param fraglen Mean RNA fragment length. Sequences will be read off the end(s) of these 
#' fragments.
#' @param fragsd Standard deviation of fragment lengths. 
#' @param readlen Read length
#' @param error_rate=0.005 Sequencing error rate. Must be between 0 and 1. A uniform error model is 
#' assumed. 
#' @param paired If \code{TRUE}, paired-end reads are simulated; else single-end reads are 
#' simulated.
#' @param seed Optional seed to set before simulating reads, for reproducibility.
#' @export
#' @details Reads can either be simulated from a FASTA file of transcripts (provided with the 
#' \code{fasta} argument) or from a GTF file plus DNA sequences (provided with the \code{gtf} and 
#' \code{seqpath} arguments). Simulating from a GTF file and DNA sequences may be a bit slower: 
#' it took about 6 minutes to parse the GTF/sequence files for chromosomes 1-22, X, and Y in hg19.
#' @examples \dontrun{
#' fastapath = system.file("data", "chr22.fa", package="polyester")
#' numtx = count_transcripts(fastapath)
#' readmat = matrix(20, ncol=10, nrow=numtx)
#' readmat[1:30, 1:5] = 40
#' 
#' simulate_experiment_countmat(fastapath, readmat, outdir="./data/", seed=5)
#'}
simulate_experiment_countmat = function(fasta=NULL, gtf=NULL, seqpath=NULL, readmat, outdir="", 
    fraglen=250, fragsd=25, readlen=100, error_rate=0.005, paired=TRUE, seed=NULL){

    if(!is.null(seed)) set.seed(seed)
    
    if(!is.null(fasta) & is.null(gtf) & is.null(seqpath)){
        transcripts = readDNAStringSet(fasta)
    }else if(is.null(fasta) & !is.null(gtf) & !is.null(seqpath)){
        message('parsing gtf and sequences...')
        transcripts = seq_gtf(gtf, seqpath, ...)
        message('done parsing')
    }else{
        stop('must provide either fasta or both gtf and seqpath')
    }

    system(paste("mkdir -p", outdir))

    for(i in 1:ncol(readmat)){
        tObj = rep(transcripts,times=readmat[,i])
  
        #get fragments
        tFrags = generate_fragments(tObj, fraglen=fraglen, fragsd=fragsd)

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