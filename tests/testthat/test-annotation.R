context('test annotation-related functions')

test_that('GenomicRanges is installed', {
    expect_that(suppressMessages(library(GenomicRanges)), not(throws_error()))
})

test_that('ballgown reader is not broken and data was installed properly', {
    suppressMessages(library(GenomicRanges))
    expect_that(ballgown(dataDir=system.file('extdata', package='ballgown'),
        samplePattern='sample', verbose=FALSE), not(throws_error()))
})

test_that('example annotation was installed properly', {
    expect_that(system.file('extdata', 'annot.gtf.gz', package='ballgown'), 
        not(equals("")))
})

## if these things fail, we should have gotten a more helpful message earlier:
bg = ballgown(dataDir=system.file('extdata', package='ballgown'), 
    samplePattern='sample', verbose=FALSE)
suppressMessages(library(GenomicRanges))
gtfPath = system.file('extdata', 'annot.gtf.gz', package='ballgown')

test_that('gtf read function works', {
    expect_that(gffRead(gtfPath), not(throws_error()))
})

x = gffRead(gtfPath)

test_that('gtf read function gives the right answer', {    
    expect_that(ncol(x), equals(9))
    expect_that(nrow(x), equals(13732))
    expect_that(names(x), is_identical_to(c('seqname', 'source', 'feature',
        'start', 'end', 'score', 'strand', 'frame', 'attributes')))
    expect_that(x, is_a('data.frame'))
})

test_that('splitting attribute fields works', {
    expect_that(getAttributeField(x$attributes, 'transcript_id'), 
        not(throws_error()))
    expect_that(all(is.na(getAttributeField(x$seqname, 'transcript_id'))), 
        is_true())
    transcripts = getAttributeField(x$attributes, 'transcript_id')
    expect_that(transcripts, is_a('character'))
    expect_that(length(transcripts), equals(13732))
    expect_that(transcripts[192], equals('ENST00000444520'))
    expect_that(all(is.na(getAttributeField(x$attributes, 'transcript_id', 
        attrsep=','))), 
        is_true())
})


test_that('reading gtf directly into GRanges and GRangesList works', {
    expect_that(gffReadGR(gtfPath), not(throws_error()))
    gr = gffReadGR(gtfPath)
    expect_that(length(gr), equals(13732))
    expect_that(gffReadGR(gtfPath, splitByTranscript=TRUE), not(throws_error()))
    grl = gffReadGR(gtfPath, splitByTranscript=TRUE)
    expect_that(length(grl), equals(910))
    expect_that(names(grl), is_a('character'))
    expect_that(sort(names(grl)), is_identical_to(
        sort(unique(getAttributeField(x$attributes, 'transcript_id')))))
})

grl = gffReadGR(gtfPath, splitByTranscript=TRUE)

test_that('pctOverlap function works', {
    expect_that(pctOverlap(structure(bg)$trans[[2]], grl[[369]]), 
        equals(0.7986441, tolerance=0.0001))
})

test_that('function for annotating assemblies works', {
    expect_that(annotate_assembly(assembled=structure(bg)$trans,
        grl), not(throws_error()))
    rel = annotate_assembly(assembled=structure(bg)$trans, grl)
    expect_that(nrow(rel), equals(1038))
    expect_that(ncol(rel), equals(3))
    expect_that(rel, is_a('data.frame'))
    expect_that(length(unique(rel$assembledInd)), equals(99))
    expect_that(length(unique(rel$annotatedInd)), equals(910))
    expect_that(rel$annotatedInd[6], equals(177))
    expect_that(rel$assembledInd[6], equals(3))
    expect_that(rel$percent[6], equals(0.8032030, tolerance=0.001))
})

test_that('function for labeling assembled transcripts with genes works', {
    expect_that(getGenes(gtfPath, structure(bg)$trans, UCSC=FALSE), 
        not(throws_error()))
    geneoverlaps = getGenes(gtfPath, structure(bg)$trans, UCSC=FALSE)
    expect_that(geneoverlaps, is_a('CharacterList'))
    expect_that(geneoverlaps[[2]], equals('ENSG00000229027'))
    expect_that(length(geneoverlaps[[4]]), equals(2))
    expect_that(length(geneoverlaps), equals(100))
})

test_that('contains function works', {
    expect_that(contains(structure(bg)$trans, grl), is_a('logical'))
    expect_that(sum(contains(structure(bg)$trans, grl)), equals(61))
    expect_that(length(contains(structure(bg)$trans, grl)), equals(100))
})

pdf(file=NULL) #don't print plots to screen while testing
test_that('plot of assembled/annotated transcripts works', {
    expect_that(checkAssembledTx(annotated=grl, assembled=structure(bg)$trans, 
        ind=4), 
        not(throws_error()))
    expect_that(checkAssembledTx(annotated=grl, 
        assembled=structure(bg)$trans[1:50], ind=2), 
        gives_warning())
    expect_that(checkAssembledTx(annotated=grl, assembled=structure(bg)$trans, 
        ind=4, 
        main='hello!'), not(throws_error()))
})
dev.off()





