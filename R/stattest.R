#' statistical tests for differential expression in ballgown
#'
#' Test each transcript, gene, exon, or intron in a ballgown object for 
#' differential expression, using comparisons of linear models.
#' 
#' @param gown name of an object of class \code{ballgown}
#' @param gowntable matrix or matrix-like object with \code{rownames} 
#'   representing feature IDs and columns representing samples, with expression
#'   estimates in the cells. Provide the feature name with \code{feature}. You
#'   must provide exactly one of \code{gown} or \code{gowntable}. NB: gowntable
#'   is log-transformed within \code{stattest} if \code{log} is \code{TRUE}, so
#'   provide un-logged expression values in \code{gowntable}.
#' @param pData Required if \code{gowntable} is provided: data frame giving
#'   phenotype data for the samples in the columns of \code{gowntable}. (Rows of
#'   \code{pData} correspond to columns of \code{gowntable}). If \code{gown} is
#'   used instead, it must have a non-null, valid \code{pData} slot (and the 
#'   \code{pData} argument to \code{stattest} should be left \code{NULL}).
#' @param mod object of class \code{model.matrix} representing the design matrix
#'   for the linear regression model including covariates of interest
#' @param mod0 object of class \code{model.matrix} representing the design 
#'   matrix for the linear regression model without the covariates of interest.
#' @param feature the type of genomic feature to be tested for differential 
#'   expression. If \code{gown} is used, must be one of \code{"gene"}, 
#'   \code{"transcript"}, \code{"exon"}, or \code{"intron"}. If \code{gowntable}
#'   is used, this is just used for labeling and can be whatever the rows of 
#'   \code{gowntable} represent.
#' @param meas the expression measurement to use for statistical tests.  Must be
#'   one of \code{"cov"}, \code{"FPKM"}, \code{"rcount"}, \code{"ucount"}, 
#'   \code{"mrcount"}, or \code{"mcov"}. Not all expression measurements are 
#'   available for all features. Leave as default if \code{gowntable} is 
#'   provided.
#' @param timecourse if \code{TRUE}, tests whether or not the expression 
#'   profiles of genomic features vary over time (or another continuous 
#'   covariate) in the study.  Default \code{FALSE}.  Natural splines are used
#'   to fit time profiles, so you must have more timepoints than degrees of 
#'   freedom used to fit the splines. The default df is 4.
#' @param df degrees of freedom used for modeling expression over time with 
#'   natural cubic splines.  Default 4.  Only used if \code{timecourse=TRUE}.
#' @param covariate string representing the name of the covariate of interest 
#'   for the differential expression tests.  Must correspond to the name of a 
#'   column of \code{pData(gown)}. If \code{timecourse=TRUE}, this should be the
#'   study's time variable.
#' @param adjustvars optional vector of strings representing the names of 
#'   potential confounders.  Must correspond to names of columns of 
#'   \code{pData(gown)}.
#' @param gexpr optional data frame that is the result of calling 
#'   \code{gexpr(gown))}.  (You can speed this function up by pre-creating 
#'   \code{gexpr(gown)}.)
#' @param getFC if \code{TRUE}, also return estimated fold changes (adjusted for
#'   library size and confounders) between populations. Only available for 
#'   2-group comparisons at the moment. Default \code{FALSE}.
#' @param libadjust library-size adjustment to use in linear models. By default,
#'   the adjustment is defined as the sum of the sample's log expression
#'   measurements below the 75th percentile of those measurements. To use 
#'   a different library-size adjustment, provide a numeric vector of each 
#'   sample's adjustment value. Entries of this vector correspond to samples in
#'   in rows of \code{pData}. If no library size adjustment is desired, set to
#'   FALSE.
#' @param log if \code{TRUE}, outcome variable in linear models is 
#'   log(expression+1), otherwise it's expression. Default TRUE.
#' 
#' @details At minimum, you need to provide a ballgown object or count table, 
#'   the type of feature you want to test (gene, transcript, exon, or intron), 
#'   the expression measurement you want to use (FPKM, cov, rcount, etc.), and 
#'   the covariate of interest, which must be the name of one of the columns of 
#'   the `pData` component of your ballgown object (or provided pData). This 
#'   covariate is automatically converted to a factor during model fitting in 
#'   non-timecourse experiments.
#' 
#'   By default, models are fit using \code{log2(meas + 1)} as the outcome for 
#'   each feature. To disable the log transformation, provide `log = FALSE` as 
#'   an argument to `stattest`. You can use the \code{gowntable} option if you'd
#'   like to to use a different transformation.
#' 
#'   Library size adjustment is performed by default by using the sum of the log
#'   nonzero expression measurements for each sample, up to the 75th percentile
#'   of those measurements. This adjustment can be disabled by setting 
#'   \code{libadjust=FALSE}. You can use \code{mod} and \code{mod0} to specify 
#'   alternative library size adjustments.
#' 
#'   \code{mod} and \code{mod0} are optional arguments.  If \code{mod} is 
#'   specified, you must also specify \code{mod0}.  If neither is specified, 
#'   \code{mod0} defaults to the design matrix for a model including only a 
#'   library-size adjustment, and \code{mod} defaults to the design matrix for a
#'   model including a library-size adjustment and \code{covariate}. Note that 
#'   if you supply \code{mod} and \code{mod0}, \code{covariate}, 
#'   \code{timecourse}, \code{adjustvars}, and \code{df} are ignored, so make 
#'   sure your covariate of interest and all appropriate confounder 
#'   adjustments, including library size, are specified in \code{mod} and 
#'   \code{mod0}. By default, the library-size adjustment is the sum of all
#'   counts below the 75th percentile of nonzero counts, on the log scale
#'   (log2 + 1). 
#' 
#'   Full model details are described in the supplement of 
#'   \url{http://biorxiv.org/content/early/2014/03/30/003665}.
#' 
#' @return data frame containing the columns \code{feature}, \code{id} 
#'   representing feature id, \code{pval} representing the p-value for testing 
#'   whether this feature was differentially expressed according to 
#'   \code{covariate}, and \code{qval}, the estimated false discovery rate 
#'   using this feature's signal strength as a significance cutoff. An 
#'   additional column, \code{fc}, is included if \code{getFC} is \code{TRUE}.
#' 
#' @export
#' 
#' @references \url{http://biorxiv.org/content/early/2014/03/30/003665}
#' 
#' @author Jeff Leek, Alyssa Frazee
#' @examples
#' data(bg)
#' 
#' # two-group comparison:
#' stat_results = stattest(bg, feature='transcript', meas='FPKM', 
#'   covariate='group')
#' 
#' # timecourse test:
#' pData(bg) = data.frame(pData(bg), time=rep(1:10, 2)) #dummy time covariate
#' timecourse_results = stattest(bg, feature='transcript', meas='FPKM', 
#'   covariate='time', timecourse=TRUE)
#' 
#' # timecourse test, adjusting for group:
#' group_adj_timecourse_results = stattest(bg, feature='transcript', 
#'   meas='FPKM', covariate='time', timecourse=TRUE, adjustvars='group')
#' 
#' # custom model matrices:
#' ### create example data:
#' set.seed(43)
#' sex = sample(c('M','F'), size=nrow(pData(bg)), replace=TRUE)
#' age = sample(21:52, size=nrow(pData(bg)), replace=TRUE)
#' 
#' ### create design matrices:
#' mod = model.matrix(~ sex + age + pData(bg)$group + pData(bg)$time)
#' mod0 = model.matrix(~ pData(bg)$group + pData(bg)$time)
#'
#' ### build model: 
#' adjusted_results = stattest(bg, feature='transcript', meas='FPKM', 
#'   mod0=mod0, mod=mod)

stattest = function(gown = NULL, gowntable = NULL, pData = NULL, mod = NULL, 
    mod0 = NULL, feature = c("gene", "exon", "intron", "transcript"), 
    meas = c("cov", "FPKM", "rcount", "ucount", "mrcount", "mcov"), 
    timecourse = FALSE, covariate = NULL, adjustvars = NULL, gexpr = NULL, 
    df = 4, getFC = FALSE, libadjust = NULL, log = TRUE){

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
        feature = match.arg(feature)
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
        if((feature == "intron") & 
            !(meas %in% c("rcount", "ucount", "mrcount"))){
            stop("introns only have rcount, ucount, and mrcount measurements")
        }
        pData = pData(gown)
        if(is.null(pData) & is.null(mod)){
            stop(.makepretty('to do statistical tests, either gown must contain
                pData or you must specify models.'))
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
    if(is.null(libadjust)){
        libadjust = apply(expr, 2, function(x){
            lognz = log2(x[x!=0] + 1)
            q3 = quantile(lognz, 0.75)
            sum(lognz[lognz<q3])
        })
    }else if(!identical(libadjust, FALSE)){
            stopifnot(is.numeric(libadjust) | identical(libadjust, FALSE))
            stopifnot(length(libadjust) == nrow(pData))
    }

    if(is.null(mod) & is.null(mod0)){
        ## by default, just test whether the given covariate is important
        colind = which(names(pData) == covariate)
        if(length(colind) == 0) stop("invalid covariate name")

        ## extract the covariate
        x = pData[,colind]
                
        ## make sure there are at least 2 reps per group:
        if(any(table(x) < 2) & !timecourse){
            stop(.makepretty('There must be at least two replicates per group.
                Make sure covariate is categorical; if continuous, consider the
                timecourse option, or specify your own models with mod and
                mod0.'))
        }

        ## make sure time variable is truly continuous:
        ## (if not, continue using it just as "important")
        if(timecourse){
            n_unique_times = length(table(x))
            if(n_unique_times <= df){
                warning(paste0('Not enough timepoints (or values of', 
                    ' covariate) to fit a spline model with ', df, ' degrees ',
                    'of freedom. Statistical tests will be run treating time',
                    ' as categorical. You can also re-run the analysis with ',
                    'decreased df.'))
                timecourse = FALSE
            }
        }

        ## create model matrices
        if(!is.null(adjustvars)){
            variable_list = ""
            for(i in seq_along(adjustvars)){
                if(!adjustvars[i] %in% names(pData)){
                    stop(paste(adjustvars[i], 'is not a valid covariate'))
                }
                column_ind = which(names(pData) == adjustvars[i])
                eval(parse(text=paste0(adjustvars[i], 
                    " <- pData[,",column_ind,"]")))
                variable_list = paste(variable_list, adjustvars[i], sep="+")
            }
            if(!identical(libadjust, FALSE)){
                eval(parse(text=paste0("mod0 = model.matrix(~ libadjust", 
                    variable_list, ")")))
                if(timecourse){
                    eval(parse(text=paste0(
                        "mod = model.matrix(~ ns(x, df = ", df, ") + libadjust", 
                        variable_list, ")")))
                } else {
                    eval(parse(text=paste0(
                        "mod = model.matrix(~ as.factor(x) + libadjust", 
                        variable_list, ")")))
                }
            } else {
                variable_list = substr(variable_list, 2, nchar(variable_list)) 
                #^^strip off "+" at beginning of variable_list
                eval(parse(text=paste0("mod0 = model.matrix(~", 
                    variable_list, ")")))
                if(timecourse){
                    eval(parse(text=paste0(
                        "mod = model.matrix(~ ns(x, df = ", df, ") + ", 
                        variable_list,")")))
                } else {
                    eval(parse(text=paste0(
                        "mod = model.matrix(~ as.factor(x) + ", variable_list, 
                        ")")))
                }                
            }
        } else {
            if(!identical(libadjust, FALSE)){
                mod0 = model.matrix(~ libadjust)
            } else{
                mod0 = matrix(1, nrow=length(x), ncol=1)
            }
            if(timecourse){
                if(!identical(libadjust, FALSE)){
                    mod = model.matrix(~ ns(x, df = df) + libadjust)
                } else {
                    mod = model.matrix(~ ns(x, df = df))
                }
            } else {
                if(!identical(libadjust, FALSE)){
                    mod = model.matrix(~ as.factor(x) + libadjust)
                } else {
                    mod = model.matrix(~ as.factor(x))
                }
            }
        }
    } else {
        ## test whether custom models are, in fact, nested.
        stopifnot(class(mod0) == 'matrix' & class(mod) == 'matrix')
        if(!('assign' %in% names(attributes(mod0))) | !('assign' %in% 
                names(attributes(mod)))){
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
                warning(.makepretty('log is FALSE, so estimated difference (not
                    fold change) is reported.'))
                estFC = lmodels$coefficients[,2]
            }
            presults = f.pvalue(y, mod, mod0)
            results = data.frame(feature=rep(feature, nrow(expr)), 
                id=rownames(expr), fc=estFC, 
                pval=presults, qval=p.adjust(presults, "fdr"))
            if(!log){
                names(results)[3] = 'difference'
            }
            rownames(results) = NULL
            return(results)
        }
    }

    presults = f.pvalue(y, mod, mod0)
    results = data.frame(feature=rep(feature, nrow(expr)), 
        id=rownames(expr), pval=presults, qval=p.adjust(presults, "fdr")) 
    rownames(results) = NULL
    return(results)
}


