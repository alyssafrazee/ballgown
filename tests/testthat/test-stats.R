context('test differential expression functions')

test_that('ballgown reader is not broken and data was installed properly', {
    expect_that(ballgown(dataDir=system.file('extdata', package='ballgown'),
        samplePattern='sample', verbose=FALSE), not(throws_error()))
})

bg = ballgown(dataDir=system.file('extdata', package='ballgown'), 
    samplePattern='sample', verbose=FALSE)
pData(bg) = data.frame(id=sampleNames(bg), group=rep(c(1,0), each=10))

test_that('default stat models work', {
    expect_that(stattest(bg, feature='transcript', meas='FPKM', 
        covariate='group'), not(throws_error()))
    expect_that(stattest(bg, feature='transcript', meas='cov', 
        covariate='group'), not(throws_error()))
    expect_that(stattest(bg, feature='transcript', meas='ucount', 
        covariate='group'), throws_error())
    statres = stattest(bg, feature='transcript', meas='FPKM', covariate='group')
    expect_that(statres, is_a('data.frame'))
    expect_that(nrow(statres), equals(100))
    expect_that(ncol(statres), equals(4))
    expect_that(statres$pval[6], equals(0.2731, tolerance=0.001))

    expect_that(stattest(bg, feature='gene', meas='FPKM', covariate='group'), 
        not(throws_error()))
    expect_that(stattest(bg, feature='gene', meas='cov', covariate='group'), 
        throws_error())
    generes = stattest(bg, feature='gene', meas='FPKM', covariate='group')
    expect_that(generes, is_a('data.frame'))
    expect_that(nrow(generes), equals(93))
    expect_that(ncol(generes), equals(4))
    expect_that(generes$pval[4], equals(0.124, tolerance=0.001))

    expect_that(stattest(bg, feature='exon', meas='cov', covariate='group'), 
        not(throws_error()))
    expect_that(stattest(bg, feature='exon', meas='FPKM', covariate='group'), 
        throws_error())
    exonres = stattest(bg, feature='exon', meas='cov', covariate='group')
    expect_that(exonres, is_a('data.frame'))
    expect_that(nrow(exonres), equals(633))
    expect_that(ncol(exonres), equals(4))
    expect_that(exonres$pval[3], equals(0.183, tolerance=0.001))

    expect_that(stattest(bg, feature='intron', meas='mrcount', 
        covariate='group'), not(throws_error()))
    expect_that(stattest(bg, feature='intron', meas='cov', covariate='group'), 
        throws_error())

    # check other errors:
    expect_that(stattest(bg, feature='junction', meas='mrcount', 
        covariate='group'), throws_error())
    expect_that(stattest(bg, feature='transcript', meas='FPKM', 
        covariate=group), throws_error())
    expect_that(stattest(bg, feature='transcript', meas='FPKM', 
        covariate='nonsense!'), throws_error())
    bgnull = bg
    pData(bgnull) = NULL
    expect_that(stattest(bgnull, feature='transcript', meas='FPKM', 
        covariate='group'), throws_error())
    rm(bgnull)
})

test_that('log, library size, and getFC options work', {
    expect_that(stattest(bg, feature='transcript', meas='FPKM', 
        covariate='group', log=FALSE), not(throws_error()))
    expect_that(stattest(bg, feature='transcript', meas='FPKM', 
        covariate='group', libadjust=FALSE), not(throws_error()))
    expect_that(stattest(bg, feature='transcript', meas='FPKM', 
        covariate='group', getFC=TRUE), not(throws_error()))
    expect_that(stattest(bg, feature='transcript', meas='FPKM', 
        covariate='group', log=FALSE, getFC=TRUE), gives_warning())
})

cpd = pData(bg)
set.seed(251)
pData(bg) = data.frame(cpd, c1=runif(20,1,5), c2=rnorm(20,0,5))

test_that('user-defined model matrices work', {
    mod0 = model.matrix(~pData(bg)$group + pData(bg)$c1)
    mod = model.matrix(~pData(bg)$group + pData(bg)$c1 + pData(bg)$c2)
    expect_that(stattest(bg, mod=mod, mod0=mod0, feature='transcript', 
        meas='FPKM'), not(throws_error()))
    customres = stattest(bg, mod=mod, mod0=mod0, feature='transcript', 
        meas='FPKM')
    expect_that(customres, is_a('data.frame'))
    expect_that(nrow(customres), equals(100))
    expect_that(ncol(customres), equals(4))
    expect_that(customres$pval[6], equals(0.27569, tolerance=0.0001))
    expect_that(stattest(bg, mod=mod, mod0=mod0, feature='transcript', 
        meas='FPKM', getFC=TRUE), gives_warning())
    expect_that(stattest(bg, mod=mod0, mod0=mod, feature='transcript', 
        meas='FPKM'), throws_error())
    m0 = model.matrix(~pData(bg)$c1 + pData(bg)$c2)
    m = model.matrix(~pData(bg)$c2 + pData(bg)$group + rnorm(20)) #non-nested
    expect_that(stattest(bg, mod=m, mod0=m0, feature='transcript', meas='FPKM'),
        throws_error())
})

test_that('adjustment for covariates works', {
    expect_that(stattest(bg, feature='transcript', meas='FPKM', 
        covariate='group', adjustvars='c1'), not(throws_error()))
    adjres = stattest(bg, feature='transcript', meas='FPKM', covariate='group', 
        adjustvars='c1')
    expect_that(adjres, is_a('data.frame'))
    expect_that(ncol(adjres), equals(4))
    expect_that(nrow(adjres), equals(100))
    expect_that(adjres$pval[1], equals(0.0218, tolerance=0.001))
    expect_that(stattest(bg, feature='transcript', meas='FPKM', 
        covariate='group', adjustvars=c('c1', 'c2')), not(throws_error()))
    expect_that(stattest(bg, feature='transcript', meas='FPKM', 
        covariate='group', adjustvars=c('c1', 'c2', 'nonsense')), 
        throws_error())
    expect_that(stattest(bg, feature='transcript', meas='FPKM', 
        covariate='group', adjustvars='group'), throws_error())
})


test_that('timecourse option works', {
    expect_that(stattest(bg, feature='transcript', meas='FPKM', covariate='c1', 
        timecourse=TRUE), not(throws_error()))
    tres = stattest(bg, feature='transcript', meas='FPKM', covariate='c1', 
        timecourse=TRUE)
    expect_that(tres, is_a('data.frame'))
    expect_that(nrow(tres), equals(100))
    expect_that(ncol(tres), equals(4))
    expect_that(tres$pval[3], equals(0.569, tolerance=0.001))

    expect_that(stattest(bg, feature='transcript', meas='FPKM', covariate='c1', 
        adjustvars='group', timecourse=TRUE), not(throws_error()))
    expect_that(stattest(bg, feature='transcript', meas='FPKM', covariate='c1'),
        throws_error())
    expect_that(stattest(bg, feature='transcript', meas='FPKM', covariate='c1',
        timecourse=TRUE, df=3), not(throws_error()))
})

test_that('stattest_table functionality works', {
    tab = texpr(bg)
    orig_res = stattest(bg, feature='transcript', meas='FPKM', 
        covariate='group')
    tab_res = stattest(gowntable=tab, pData=pData(bg), feature='transcript', 
        covariate='group')
    expect_that(tab_res, is_identical_to(orig_res))
})