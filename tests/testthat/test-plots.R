context('test plotting functions')

bg = ballgown(dataDir=system.file('extdata', package='ballgown'), 
    samplePattern='sample', verbose=FALSE)
pData(bg) = data.frame(id=sampleNames(bg), group=rep(c(1,0), each=10))

pdf(file=NULL) # don't print plots to screen while testing

test_that('plotTranscripts does not throw errors', {
    # multi-transcript gene:
    expect_that(plotTranscripts('XLOC_000454', bg), not(throws_error()))
    expect_that(plotTranscripts('XLOC_000454', bg, samples='sample06'), 
        not(throws_error()))
    expect_that(plotTranscripts('XLOC_000454', bg, 
        samples=c('sample06', 'sample08')), not(throws_error()))
    expect_that(plotTranscripts('XLOC_000454', bg, colorby='exon', meas='cov'),
        not(throws_error()))
    expect_that(plotTranscripts('XLOC_000454', bg, colorby='exon'), 
        throws_error())
    expect_that(plotTranscripts('XLOC_000454', bg, colorby='transcript',
        customCol=c('red', 'blue', 'green'), legend=FALSE), 
    not(throws_error()))

    # single-transcript gene:
    expect_that(plotTranscripts('XLOC_000010', bg), not(throws_error()))
})

test_that('plotMeans does not throw errors', {
    expect_that(plotMeans('XLOC_000454', bg, groupvar='group'), 
        not(throws_error()))
    expect_that(plotMeans('XLOC_000454', bg, groupvar='group'), 
        not(throws_error()))
    expect_that(plotMeans('XLOC_000454', bg, groupvar='group', groupname=0), 
        not(throws_error()))
    expect_that(plotMeans('XLOC_000454', bg, groupvar='group', overall=TRUE), 
        gives_warning())
    expect_that(plotMeans('XLOC_000454', bg, groupvar='group', overall=FALSE), 
        not(throws_error()))
    expect_that(plotMeans('XLOC_000454', bg, groupvar='group', colorby='exon', 
        meas='cov'), not(throws_error()))
})

test_that('plotLatentTranscripts does not throw errors', {
    expect_that(plotLatentTranscripts('XLOC_000454', bg, method='kmeans', k=2), 
        not(throws_error()))
    expect_that(plotLatentTranscripts('XLOC_000454', bg, method='kmeans', k=8), 
        throws_error())
    expect_that(plotLatentTranscripts('XLOC_000454', bg, method='kmeans', 
        choosek='var90'), not(throws_error()))
    expect_that(plotLatentTranscripts('XLOC_000454', bg, method='hclust', k=2), 
        not(throws_error()))
}) #this function is mostly clusterTranscripts anyway

dev.off()