#' simulate RNA-seq experiment using negative binomial model
#'
#' create FASTA files containing RNA-seq reads simulated from provided transcripts, with optional differential expression between two groups
#' @param fasta path to FASTA file containing transcripts from which to simulate reads 
#' @param  num_reps How many biological replicates should be in each group? If \code{num_reps} 
#' is an integer, \code{num_reps} replicates will be simulated in each group. 
#' Otherwise, \code{num_reps} can be a length-2 vector, where \code{num_reps[1]} 
#' and \code{num_reps[2]} replicates will be simulated in each of the two groups.
#' @param  fraglen Mean RNA fragment length. Sequences will be read off the end(s) of these fragments.
#' @param  fragsd Standard deviation of fragment lengths. 
#' @param  readlen Read length.
#' @param  error_rate Sequencing error rate. Must be between 0 and 1. A uniform error model is assumed. 
#' @param  paired If \code{TRUE}, paired-end reads are simulated; else single-end reads are simulated.
#' @param  reads_per_transcript baseline mean number of reads to simulate 
#' from each transcript. Can be an integer, in which case this many reads
#' are simulated from each transcript, or an integer vector whose length
#' matches the number of transcripts in \code{fasta}. 
#' @param  fold_changes Vector of multiplicative fold changes between groups,
#' one entry per transcript in \code{fasta}. A fold change > 1 means the 
#' transcript is overexpressed in the first \code{num_reps} (or \code{num_reps[1]})
#' samples. Fold change < 1 means transcript is overexpressed in the last
#' \code{num_reps} (or \code{num_reps[2]}) samples. The change is in the mean
#' number of reads generated from the transcript, between groups.
#' @param  dispersion_param the negative binomial \code{size} distribution (see \code{\link{NegBinomial}})
#' of the number of reads drawn per transcript. If left blank, defaults to 
#' \code{reads_per_transcript / 3}. Negative binomial variance is mean + mean^2 / size.
#' @param  outdir character, path to folder where simulated reads should be written. Should end with "/" if specified. If unspecified, reads are written to current working directory.
#' @export
#' @examples \dontrun{
#' ## simulate a few reads from chromosome 22
#' 
#' fastapath = system.file("data", "chr22.fa", package="polyester")
#' numtx = count_transcripts(fastapath)
#' set.seed(4)
#' fold_changes = sample(c(0.5, 1, 2), size=numtx, probs=c(0.05, 0.9, 0.05), replace=TRUE)
#' 
#' simulate_experiment(fastapath, reads_per_transcript = 10, fold_changes=fold_changes, outdir="./simdata/")
#'
#'}
simulate_experiment = function(fasta, num_reps=10, fraglen=250, fragsd=25, 
    readlen=100, error_rate=0.005, paired=TRUE, reads_per_transcript=300, 
    fold_changes, dispersion_param=NULL, outdir=""){

    require(Biostrings)
    transcripts = readDNAStringSet(fasta)
    L = width(transcripts)
        
    ## get number of reps per group
    stopifnot(length(num_reps)==1 | length(num_reps)==2)
    if(length(num_reps)==2){
        n1 = num_reps[1]
        n2 = num_reps[2]
    }else{
        n1 = n2 = num_reps
    }

    ## get baseline means for each group from fold changes:
    if(exp(mean(log(fold_changes))) != 1){
        warning('provided fold changes mean that one group will have
            more overall expression than the other, so some of
            the DE signal might be lost due to library size adjustment.')
    }
    basemeans1 = ifelse(fold_changes > 1, 
        reads_per_transcript*fold_changes, reads_per_transcript)
    basemeans1 = round(basemeans1)
    basemeans2 = ifelse(fold_changes < 1,
        reads_per_transcript/fold_changes, reads_per_transcript)
    basemeans2 = round(basemeans2)

    if(is.null(dispersion_param)){
        dispersion_param = reads_per_transcript/3
    }

    if(length(dispersion_param) == 1){
        nbdp1 = nbdp2 = dispersion_param
    }else{
        nbdp1 = dispersion_param[1:n1]
        nbdp2 = dispersion_param[(1:n2)+n1] 
    }

    numreadsList = vector("list", n1+n2)
    for(i in 1:n1){
        numreadsList[[i]] = NB(basemeans1, nbdp1)
    }
    for(i in (1:n2)+n1){
        numreadsList[[i]] = NB(basemeans2, nbdp2)
    }


    ## simulate reads for each sample:
    #######################################
    for(i in 1:(n1+n2)){
        
        tObj = rep(transcripts, times=numreadsList[[i]])
        
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
            fname=paste0(outdir, 'sample_', sprintf('%02d', i)), paired=paired)
    }
}
