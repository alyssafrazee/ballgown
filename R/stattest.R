#' statistical tests for differential expression in ballgown
#'
#' Test each transcript, gene, exon, or intron in a ballgown object for differential expression, 
#' using comparisons of linear models.
#' 
#' @param gown name of an object of class \code{ballgown}
#' @param gowntable matrix or matrix-like object with \code{rownames} representing feature IDs and 
#'   columns representing samples, with expression estimates in the cells. Provide the feature name 
#'   with \code{feature}. You must provide exactly one of \code{gown} or \code{gowntable}.
#' @param pData Required if \code{gowntable} is provided: data frame giving phenotype data for the
#'    samples in the columns of \code{gowntable}. (Rows of \code{pData} correspond to columns of
#'    \code{gowntable}). If \code{gown} is used instead, it must have a non-null, valid \code{pData}
#'    slot (and the \code{pData} argument to \code{stattest} should be left \code{NULL}).
#' @param mod object of class \code{model.matrix} representing the design matrix for the linear 
#'   regression model including covariates of interest
#' @param mod0 object of class \code{model.matrix} representing the design matrix for the linear 
#'   regression model without the covariates of interest.
#' @param feature the type of genomic feature to be tested for differential expression. Must be one 
#'   of \code{"gene"}, \code{"transcript"}, \code{"exon"}, or \code{"intron"}.
#' @param meas the expression measurement to use for statistical tests.  Must be one of 
#'   \code{"cov"}, \code{"FPKM"}, \code{"rcount"}, \code{"ucount"}, \code{"mrcount"}, or 
#'   \code{"mcov"}. Not all expression measurements are available for all features. Leave as default
#'    if \code{gowntable} is provided.
#' @param timecourse if \code{TRUE}, tests whether or not the expression profiles of genomic 
#'   features vary over time (or another continuous covariate) in the study.  Default \code{FALSE}.
#' @param covariate string representing the name of the covariate of interest for the differential 
#'   expression tests.  Must correspond to the name of a column of \code{pData(gown)}. If 
#'   \code{timecourse=TRUE}, this should be the study's time variable.
#' @param adjustvars optional vector of strings representing the names of potential confounders.  
#'   Must correspond to names of columns of \code{pData(gown)}.
#' @param gexpr optional data frame that is the result of calling \code{gexpr(gown))}.  (You can 
#'   speed this function up by pre-creating \code{gexpr(gown)}.)
#' @param df degrees of freedom used for modeling expression over time with natural cubic splines.  
#'   Default 4.  Only used if \code{timecourse=TRUE}.
#' @param getFC if \code{TRUE}, also return estimated fold changes (adjusted for library size and 
#'   confounders) between populations. Only available for 2-group comparisons at the moment. Default 
#'   \code{FALSE}.
#' @param libadjust if \code{TRUE} (default), include a library-size adjustment as a confounder in 
#'   the fitted models. The adjustment is currently defined as the sum of the sample's counts below 
#'   that sample's 75th percentile. You can change the library size adjustment by setting
#'   \code{libadjust} to \code{FALSE} and providing \code{mod} and \code{mod0} with the desired
#'   library size adjustments.
#' @param log if \code{TRUE}, outcome variable in linear models is log(expression+1), otherwise it's 
#'   expression. Default TRUE.
#' 
#' @details \code{mod} and \code{mod0} are optional arguments.  If \code{mod} is specified, you must
#'   also specify \code{mod0}.  If neither is specified, \code{mod0} defaults to the design matrix 
#'   for a model including only a library-size adjustment, and \code{mod} defaults to the design 
#'   matrix for a model including a library-size adjustment and \code{covariate}. Note that if you 
#'   supply \code{mod} and \code{mod0}, \code{covariate}, \code{timecourse}, \code{adjustvars}, and 
#'   \code{df} are ignored, so make sure your covariate of interest and all appropriate confounder 
#'   adjustments, including library size, are specified in \code{mod} and \code{mod0}.
#' 
#' @return data frame containing the columns \code{feature}, \code{id} representing feature id, 
#'   \code{pval} representing the p-value for testing whether this feature was differentially 
#'   expressed according to \code{covariate}, and \code{qval}, the estimated false discovery rate 
#'   using this feature's signal strength as a significance cutoff. An additional column, \code{fc}, 
#'   is included if \code{getFC} is \code{TRUE}.

#' @export

#' @author Jeff Leek, Alyssa Frazee

stattest = function(gown = NULL, gowntable = NULL, pData = NULL, mod = NULL, mod0 = NULL, 
    feature = c("gene", "exon", "intron", "transcript"), 
    meas = c("cov", "FPKM", "rcount", "ucount", "mrcount", "mcov"), timecourse = FALSE, 
    covariate = NULL, adjustvars = NULL, gexpr = NULL, df = 4, getFC = FALSE, libadjust = TRUE, 
    log = TRUE){

    feature = match.arg(feature)

    if(!xor(is.null(gown), is.null(gowntable))){
        stop('must provide exactly one of gown and gowntable')
    }
    if(xor(is.null(mod), is.null(mod0))){
        stop("please provide both null and full models, or use the defaults")
    }

    if(!is.null(gowntable)){
        expr = as.matrix(gowntable)
        stopifnot(!is.null(pData))
        stopifnot(nrow(pData) == ncol(expr))
    } else {
        meas = match.arg(meas)
        if(feature == "transcript" & !(meas %in% c("cov", "FPKM"))){
            stop("transcripts only have cov and FPKM measurements")
        }
        if(feature == "gene" & meas != "FPKM"){
            stop("gene tests can only be done on FPKM measurements")
        }
        if((feature == "exon") & meas == "FPKM"){
            stop("exons do not have FPKM measurements")
        }
        if((feature == "intron") & !(meas %in% c("rcount", "ucount", "mrcount"))){
            stop("introns only have rcount, ucount, and mrcount measurements")
        }
        pData = pData(gown)
        if(is.null(pData) & is.null(mod)){
            stop(.makepretty('to do statistical tests, either gown must contain pData or you must
                specify models.'))
        }
            ## extract the right expression measurements
        if(feature == "gene"){
            if(is.null(gexpr)){
                expr = gexpr(gown)
            }else{
                expr = gexpr
            }
        }
        if(feature == "exon") expr = eexpr(gown, meas)
        if(feature == "intron") expr = iexpr(gown, meas)
        if(feature == "transcript") expr = texpr(gown, meas)
        expr = as.matrix(expr)
    }
    
    n = ncol(expr)

    ## library size adjustment
    if(libadjust){
        lib_adj = apply(expr, 2, function(x){
            q3 = quantile(x, 0.75)
            sum(x[x<q3])
        })
    }

    if(is.null(mod) & is.null(mod0)){
        ## by default, just test whether the given covariate is important
        colind = which(names(pData) == covariate)
        if(length(colind) == 0) stop("invalid covariate name")

        ## extract the covariate
        x = pData[,colind]
        if(length(unique(as.factor(x))) == length(x) & !timecourse){
            stop(.makepretty('Default models only support continuous covariates of interest when
                timecourse is TRUE. Consider the timecourse option, or you can specify your own
                models with mod and mod0.'))
        }
        
        ## create model matrices
        if(!is.null(adjustvars)){
            variable_list = ""
            for(i in seq_along(adjustvars)){
                if(!adjustvars[i] %in% names(pData)){
                    stop(paste(adjustvars[i], 'is not a valid covariate'))
                }
                column_ind = which(names(pData) == adjustvars[i])
                eval(parse(text=paste0(adjustvars[i], " <- pData[,",column_ind,"]")))
                variable_list = paste(variable_list, adjustvars[i], sep="+")
            }
            if(libadjust){
                eval(parse(text=paste0("mod0 = model.matrix(~ lib_adj", variable_list, ")")))
                if(timecourse){
                    eval(parse(text=paste0("mod = model.matrix(~ ns(x, df = ", df, ") + lib_adj", 
                        variable_list, ")")))
                } else {
                    eval(parse(text=paste0("mod = model.matrix(~ as.factor(x) + lib_adj", 
                        variable_list, ")")))
                }
            } else {
                variable_list = substr(variable_list, 2, nchar(variable_list)) 
                #^^strip off "+" at beginning of variable_list
                eval(parse(text=paste0("mod0 = model.matrix(~", variable_list, ")")))
                if(timecourse){
                    eval(parse(text=paste0("mod = model.matrix(~ ns(x, df = ", df, ") + ", 
                        variable_list,")")))
                } else {
                    eval(parse(text=paste0("mod = model.matrix(~ as.factor(x) + ", variable_list, 
                        ")")))
                }                
            }
        } else {
            if(libadjust){
                mod0 = model.matrix(~ lib_adj)
            } else{
                mod0 = matrix(1, nrow=length(x), ncol=1)
            }
            if(timecourse){
                if(libadjust){
                    mod = model.matrix(~ ns(x, df = df) + lib_adj)
                } else {
                    mod = model.matrix(~ ns(x, df = df))
                }
            } else {
                if(libadjust){
                    mod = model.matrix(~ as.factor(x) + lib_adj)
                } else {
                    mod = model.matrix(~ as.factor(x))
                }
            }
        }
    } else {
        ## test whether custom models are, in fact, nested.
        stopifnot(class(mod0) == 'matrix' & class(mod) == 'matrix')
        if(!('assign' %in% names(attributes(mod0))) | !('assign' %in% names(attributes(mod)))){
            stop('mod and mod0 must both be model.matrix objects')
        }
        stopifnot(ncol(mod) > ncol(mod0))
        for(i in 1:ncol(mod0)){
            inMod = any(apply(mod, 2, function(x) identical(mod0[,i], x)))
            if(!inMod){
                stop('mod0 is not nested in mod')
            }
        }
    }

    if(log){
        y = log2(expr+1)
    }else{
        y = expr
    }

    if(getFC){
        two = try(length(table(x)), silent=TRUE)
        if(class(two) == 'try-error'){
            warning('fold changes not available for custom models')
        } else if(two != 2){
            warning('fold changes only available for 2-group comparisons')
        }else{
            lmodels = lmFit(y, design=mod)
            if(log){
                estFC = 2^(lmodels$coefficients[,2])
            }else{
                warning(.makepretty('log is FALSE, so estimated difference (not fold change) is
                    reported.'))
                estFC = lmodels$coefficients[,2]
            }
            presults = f.pvalue(y, mod, mod0)
            results = data.frame(feature=rep(feature, nrow(expr)), id=rownames(expr), fc=estFC, 
                pval=presults, qval=p.adjust(presults, "fdr"))
            if(!log){
                names(results)[3] = 'difference'
            }
            rownames(results) = NULL
            return(results)
        }
    }

    presults = f.pvalue(y, mod, mod0)
    results = data.frame(feature=rep(feature, nrow(expr)), id=rownames(expr), pval=presults, 
        qval=p.adjust(presults, "fdr")) 
    rownames(results) = NULL
    return(results)
}


