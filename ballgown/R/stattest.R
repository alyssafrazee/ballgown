#' statistical tests for differential expression in ballgown
#'
#' @param gown name of an object of class \code{ballgown}
#' @param mod object of class \code{model.matrix} representing the design matrix for the linear regression model including covariates of interest
#' @param mod0 object of class \code{model.matrix} representing the design matrix for the linear regression model without the covariates of interest.
#' @param feature the type of genomic feature to be tested for differential expression. Must be one of \code{"gene"}, \code{"transcript"}, \code{"exon"}, or \code{"intron"}.
#' @param meas the expression measurement to use for statistical tests.  Must be one of \code{"cov"}, \code{"FPKM"}, \code{"rcount"}, \code{"ucount"}, \code{"mrcount"}, or \code{"mcov"}. Not all expression measurements are available for all features.
#' @param timecourse if \code{TRUE}, tests whether or not the expression profiles of genomic features vary over time in the study.  Default \code{FALSE}.
#' @param covariate string representing the name of the covariate of interest for the differential expression tests.  Must correspond to the name of a column of \code{pData(gown)}. If \code{timecourse=TRUE}, this should be the study's time variable.
#' @param adjustvars optional vector of strings representing the names of potential confounders.  Must correspond to names of columns of \code{pData(gown)}.
#' @param gexpr optional data frame that is the result of calling \code{gexp(gown))}.  Used to speed this up by pre-creating \code{gexp(gown)}.
#' @param df degrees of freedom used for modeling expression over time with natural cubic splines.  Default 4.  Only used if \code{timecourse=TRUE}.
#' @details \code{mod} and \code{mod0} are optional arguments.  If \code{mod} is specified, you must also specify \code{mod0}.  If neither are specified, \code{mod0} defaults to the design matrix for a model including only a library-size adjustment, and \code{mod} defaults to the design matrix for a model including a library-size adjustment and \code{covariate}.
#' @return data frame containing the columns \code{feature}, \code{id} representing feature id, \code{pval} representing the p-value for testing whether this feature was differentially expressed according to \code{covariate}, and \code{qval}, the estimated false discovery rate using this feature's signal strength as a significance cutoff.
#' @export
#' @author Alyssa Frazee, Jeff Leek

stattest = function(gown, mod = NULL, mod0 = NULL, 
    feature = c("gene", "exon", "intron", "transcript"),
    meas = c("cov", "FPKM", "rcount", "ucount", "mrcount", "mcov"),
    timecourse = FALSE,
    covariate = NULL,
    adjustvars = NULL,
    gexpr = NULL,
    df = 4){

    feature = match.arg(feature)
    meas = match.arg(meas)

    ## check input
    if(feature == "transcript" & 
        !(meas %in% c("cov", "FPKM"))){
        stop("transcripts only have cov and FPKM measurements")
    }
    if(feature == "gene" & meas != "FPKM"){
        stop("gene tests can only be done on FPKM measurements")
    }
    if((feature == "exon" | feature == "intron") & meas == "FPKM"){
        stop("exons and introns do not have FPKM measurements")
    }
    if(xor(is.null(mod), is.null(mod0))){
        stop("please provide both null and full models, or use the defaults")
    }

    ## extract the right expression measurements
    if(feature == "gene"){
        if(is.null(gexpr)){
            gnames = indexes(gown)$t2g$g_id
            inds_by_gene = split(seq(along=gnames), gnames)
            tmeas = texpr(gown, "FPKM")
            gid_by_exon = lapply(1:nrow(texpr(gown)), function(i){rep(texpr(gown)$gene_id[i], texpr(gown)$num_exons[i])})
            ulstruct = unlist(structure(gown)$trans)
            glist = split(ulstruct, unlist(gid_by_exon))
            glengths = sapply(width(reduce(glist)), sum)
            tlengths = sapply(width(structure(gown)$trans), sum)
            tfrags = lapply(1:nrow(tmeas), function(i){
                (tlengths[i]/1000) * tmeas[i,]
            }) ## 7 minutes, too slow still but manageable
            tfrags = matrix(unlist(tfrags, use.names=FALSE), nrow=length(tfrags), byrow=TRUE)
            expr = t(sapply(1:length(inds_by_gene), function(i){colSums(tfrags[inds_by_gene[[i]],,drop=FALSE]) / glengths[i]}))
            rownames(expr) = names(inds_by_gene)
        }else{
            expr = gexpr
        }
    }
    if(feature == "exon") expr = eexpr(gown, meas)
    if(feature == "intron") expr = iexpr(gown, meas)
    if(feature == "transcript") expr = texpr(gown, meas)

    expr = as.matrix(expr)
    n = ncol(expr)

    ## library size adjustment
    lib_adj = apply(expr, 2, function(x){
        q3 = quantile(x, 0.75)
        sum(x[x<q3])
    })


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
                eval(parse(text=paste0(adjustvars[i], " <- pData(gown)[,",column_ind,"]")))
                variable_list = paste(variable_list, adjustvars[i], sep="+")
            }
            eval(parse(text=paste0("mod0 = model.matrix(~ lib_adj",variable_list,")")))
            if(timecourse){
                eval(parse(text=paste0("mod = model.matrix(~ ns(x, df = ",df,") + lib_adj",variable_list,")")))
            } else {
                eval(parse(text=paste0("mod = model.matrix(~ as.factor(x) + lib_adj",variable_list,")")))
            }
        } else {
            mod0 = model.matrix(~ lib_adj)
            if(timecourse){
                mod = model.matrix(~ ns(x, df = df) + lib_adj)
            } else {
                mod = model.matrix(~ as.factor(x) + lib_adj)
            }
        }
    }

    results = f.pvalue(log2(expr + 1), mod, mod0)
    return(data.frame(feature=rep(feature, nrow(expr)), id=rownames(expr), pval = results, qval=p.adjust(results, "fdr")))
}

