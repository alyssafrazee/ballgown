context('test ballgown constructor + slots')

## need to load object first
## this will test constructor, pData, pData<-, and sampleNames for errors...
bg = ballgown(dataDir=system.file('extdata', package='ballgown'), 
    samplePattern='sample', verbose=FALSE)
pData(bg) = data.frame(id=sampleNames(bg), group=rep(c(1,0), each=10))
suppressMessages(library(GenomicRanges))

test_that('all slots are present and have correct class', {
    expect_that(bg@expr, is_a('list'))
    expect_that(bg@structure, is_a('list'))
    expect_that(bg@indexes, is_a('list'))

    expect_that(bg@expr$exon, is_a('data.frame'))
    expect_that(eexpr(bg), is_a('matrix'))
    expect_that(bg@expr$trans, is_a('data.frame'))
    expect_that(texpr(bg), is_a('matrix'))
    expect_that(bg@expr$intron, is_a('data.frame'))
    expect_that(iexpr(bg), is_a('matrix'))

    expect_that(bg@structure$exon, is_a('GRanges'))
    expect_that(bg@structure$intron, is_a('GRanges'))
    expect_that(bg@structure$trans, is_a('GRangesList'))

    expect_that(bg@indexes$e2t, is_a('data.frame'))
    expect_that(bg@indexes$i2t, is_a('data.frame'))
    expect_that(bg@indexes$t2g, is_a('data.frame'))
    expect_that(bg@indexes$pData, is_a('data.frame'))
    expect_that(bg@indexes$bamfiles, is_identical_to(NULL))

    expect_that(sampleNames(bg), is_a('character'))
})

test_that('expression data is read accurately', {
    expect_that(dim(texpr(bg)), is_equivalent_to(c(100, 20)))
    expect_that(ncol(texpr(bg, 'all')), equals(50))
    expect_that(dim(eexpr(bg)), is_equivalent_to(c(633, 20)))
    expect_that(ncol(eexpr(bg, 'all')), equals(145))
    expect_that(dim(iexpr(bg)), is_equivalent_to(c(536, 20)))
    expect_that(ncol(iexpr(bg, 'all')), equals(65))
    expect_that(texpr(bg)[29, 11], equals(171.12))
    expect_that(eexpr(bg)[9, 19], equals(21))
    expect_that(iexpr(bg, 'ucount')[89, 2], equals(8))
})

test_that('indexes are read accurately', {
    expect_that(dim(indexes(bg)$t2g), is_equivalent_to(c(100, 2)))
    expect_that(indexes(bg)$t2g$t_id[19], equals(258))
    expect_that(dim(indexes(bg)$e2t), is_equivalent_to(c(684, 2)))
    expect_that(indexes(bg)$e2t$t_id[131], equals(377))
    expect_that(dim(indexes(bg)$i2t), is_equivalent_to(c(584, 2)))
    expect_that(indexes(bg)$i2t$i_id[512], equals(2627))
})

test_that('structure information is read accurately', {
    expect_that(length(structure(bg)$intron), equals(536))
    expect_that(start(ranges(structure(bg)$intron))[56], equals(22285664))
    expect_that(length(structure(bg)$exon), equals(633))
    expect_that(end(ranges(structure(bg)$exon))[194], equals(29747255))
    expect_that(length(structure(bg)$trans), equals(100))
    expect_that(width(ranges(structure(bg)$trans[[7]]))[3], equals(159))
})















