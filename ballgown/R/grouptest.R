# statistical test for group comparisons in ballgown
# handles 2+ groups

grouptest = function(gown, mod = NULL, mod0 = NULL, 
	feature = c("gene", "exon", "intron", "transcript"),
    expression_meas = c("cov", "FPKM", "rcount", "ucount", "mrcount", "mcov"),
	groupvar, adjustvars = NULL){

	feature = match.arg(feature)
    expression_meas = match.arg(expression_meas)

    if((feature == "gene" | feature == "transcript") & 
        !(expression_meas %in% c("cov", "FPKM"))){
        stop("genes/transcripts only have cov and FPKM measurements")
    }
    if((feature == "exon" | feature == "intron") & expression_meas == "FPKM"){
        stop("exons and introns do not have FPKM measurements")
    }

	if(feature == "gene"){
        gnames = indexes(gown)$t2g$g_id
        inds_by_gene = split(seq(along=gnames), gnames)
        expr = t(sapply(inds_by_gene, function(i){
        	colSums(texpr(gown)[i,-c(1:10),drop=FALSE])}
        ))
	}
	if(feature == "exon") expr = eexpr(gown)[,-c(1:5)]
	if(feature == "intron") expr = iexpr(gown)[,-c(1:5)]
	if(feature == "transcript") expr = texpr(gown)[,-c(1:10)]

    ## subset to only the right measurement and convert to matrix
    expr = as.matrix(expr[,grepl(expression_meas, colnames(expr))])

	n = ncol(expr)

    ## library size adjustment
	med_expr = colMedians(expr)

    if(is.null(mod) & is.null(mod0)) {
    	## just run simple test: is the group variable important?
        group_col_ind = which(names(pData(gown)) == group)
        if(length(group_col_ind) == 0) stop("invalid group variable name")

        ## create group covariate
        group = pData(gown)[,group_col_ind]
        
        ## create model matrices
        if(!is.null(adjustvars)){
            variable_list = ""
            for(i in seq_along(adjustvars)){
                column_ind = which(names(pData(gown)) == adjustvars[i])
                eval(parse(text=paste0(adjustvars[i], " <- pData(gown)[,",i,"]")))
                variable_list = paste(variable_list, adjustvars[i], sep="+")
            }
            eval(parse(text=paste0("mod0 = model.matrix(~ med_expr",variable_list,")")))
            eval(parse(text=paste0("mod = model.matrix(~ as.factor(group) + med_expr",variable_list,")")))
        } else {
            mod0 = model.matrix(~ med_expr)
            mod = model.matrix(~ as.factor(group) + med_expr)
        }
    }

    results = f.pvalue(log2(expr + 1), mod, mod0)

    return(data.frame(feature=rep(feature, nrow(expr)), id=rownames(expr), pval = results, qval = p.adjust(results, "fdr")))

}