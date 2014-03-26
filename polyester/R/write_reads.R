# write out fasta files

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
