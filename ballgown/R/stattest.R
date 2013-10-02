# statistical tests for differential expression in ballgown

stattest = function(gown, mod = NULL, mod0 = NULL, 
	feature = c("gene", "exon", "intron", "transcript"),
    meas = c("cov", "FPKM", "rcount", "ucount", "mrcount", "mcov"),
    timecourse = FALSE,
    covariate = NULL,
    adjustvars = NULL,
    df = 4){

	feature = match.arg(feature)
    expression_meas = match.arg(expression_meas)

    ## check input
    if((feature == "gene" | feature == "transcript") & 
        !(meas %in% c("cov", "FPKM"))){
        stop("genes/transcripts only have cov and FPKM measurements")
    }
    if((feature == "exon" | feature == "intron") & meas == "FPKM"){
        stop("exons and introns do not have FPKM measurements")
    }
    if(xor(is.null(mod), is.null(mod0))){
        stop("please provide both null and full models, or use the defaults")
    }

    ## extract the right expression measurements
	if(feature == "gene"){
        gnames = indexes(gown)$t2g$g_id
        inds_by_gene = split(seq(along=gnames), gnames)
        expr = t(sapply(inds_by_gene, function(i){
        	colSums(texpr(gown, meas)[i,,drop=FALSE])}
        ))
	}
	if(feature == "exon") expr = eexpr(gown, meas)
	if(feature == "intron") expr = iexpr(gown, meas)
	if(feature == "transcript") expr = texpr(gown, meas)

    expr = as.matrix(expr)
	n = ncol(expr)

    ## library size adjustment
	med_expr = colMedians(expr)

    if(is.null(mod) & is.null(mod0)){
        ## by default, just test whether the given covariate is important
        colind = which(names(pData(gown)) == covariate)
        if(length(colind) == 0) stop("invalid covariate name")

        ## extract the covariate
        x = pData(gown)[,colind]
        
        ## create model matrices
        if(!is.null(adjustvars)){
            variable_list = ""
            for(i in seq_along(adjustvars)){
                column_ind = which(names(pData(gown)) == adjustvars[i])
                eval(parse(text=paste0(adjustvars[i], " <- pData(gown)[,",i,"]")))
                variable_list = paste(variable_list, adjustvars[i], sep="+")
            }
            eval(parse(text=paste0("mod0 = model.matrix(~ med_expr",variable_list,")")))
            if(timecourse){
                eval(parse(text=paste0("mod = model.matrix(~ ns(x, df = ",df,") + med_expr",variable_list,")")))
            } else {
                eval(parse(text=paste0("mod = model.matrix(~ as.factor(x) + med_expr",variable_list,")")))
            }
        } else {
            mod0 = model.matrix(~ med_expr)
            if(timecourse){
                mod = model.matrix(~ ns(x, df = df) + med_expr)
            } else {
                mod = model.matrix(~ as.factor(x) + med_expr)
            }
        }
    }

    results = f.pvalue(log2(expr + 1), mod, mod0)
    return(data.frame(feature=rep(feature, nrow(expr)), id=rownames(expr), pval = results, qval=p.adjust(results, "fdr")))
}