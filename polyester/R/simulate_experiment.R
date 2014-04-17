# main wrapper function for the whole thing


fasta="../data/chr22.fa"; num_reps=10; 
fraglen=250; fragsd=25; readlen=100; error_rate=0.005;
paired=TRUE; reads_per_transcript=300; 
dispersion_param=NULL; outdir=""

simulate_experiment = function(fasta, num_reps=10, 
    fraglen=250, fragsd=25, readlen=100, error_rate=0.005,
    paired=TRUE, reads_per_transcript=300, fold_changes,
    dispersion_param=NULL, outdir=""){

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
        warning('provided fold changes mean that one group will have \n
            more overall expression than the other, so some of \n
            the DE signal might be lost due to library size adjustment.')
    }
    basemeans1 = ifelse(fold_changes > 1, 
        reads_per_transcript*fold_changes, reads_per_transcript)
    basemeans1 = round(basemeans1)
    basemeans2 = ifelse(fold_changes < 1,
        reads_per_transcript/fold_changes, reads_per_transcript)
    basemeans2 = round(basemeans2)

    if(is.null(dispersion_param)){
        dispersion_param = round(reads_per_transcript/3)
    }

    numreadsList = vector("list", n1+n2)
    for(i in 1:n1){
        numreadsList[[i]] = NB(1:length(transcripts), basemeans1, dispersion_param)
    }
    for(i in (1:n2)+n1){
        numreadsList[[i]] = NB(1:length(transcripts), basemeans2, dispersion_param)
    }

    ## simulate reads for each sample:
    #######################################
    for(i in 1:(n1+n2)){
        
        tObj = rep(transcripts,times=numreadsList[[i]])
        
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
