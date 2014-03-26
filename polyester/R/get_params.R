#' estimate zero-inflated negative binomial parameters from a real dataset
#'
#' @param bg ballgown object created from real RNA-seq dataset
#' @param ntranscripts number of transcripts you would like to simulate from
#' @param threshold only estimate parameters from transcripts with mean FPKM measurements larger than \code{threshold}
#' @param seed optional seed to set (for reproducibility)
#' @return list with components \code{index} (the indexes of the randomly selected transcripts from \code{bg}), \code{modfit} (the mean/dispersion model fit to the data with loess), and \code{params} (per-transcript zero-inflated negative binomial parameters \code{mu}, \code{size}, and \code{p0} for each transcript).
#' @export
#' @author Jeff Leek

get_params = function(bg,ntranscripts,threshold=3,seed=NULL){
    require(genefilter)
    if(!is.null(seed)){set.seed(seed)}
    tmeas = as.matrix(texpr(bg, "FPKM"))
    trowm = rowMeans(tmeas)
    index1 = which(trowm > threshold)
  
    tlengths = sapply(width(structure(bg)$trans[index1]), sum)
    scounts = tlengths*tmeas[index1,]/1000
    scounts = round(scounts*100)

    nsamples = dim(scounts)[2]
    scounts0 = scounts==0
    nn0 = rowSums(!scounts0)
    mu = rowSums((!scounts0)*scounts)/nn0
    s2 = rowSums((!scounts0)*(scounts - mu)^2)/(nn0-1)
    size = mu^2/(s2-mu)
    p0 = (nsamples-nn0)/nsamples
    params = data.frame(mu=mu,size=size,p0)
  
    if(ntranscripts < length(index1)){
        index2 = sample(1:length(index1),size=ntranscripts)
        index = index1[index2]
        params = params[index2,]
        scounts = scounts[index2,]
    }else{
        stop("Not enough transcripts with sufficient expression, either lower threshold or lower the number of transcripts.")
    }

    params$lsize = log2(params$size)
    params$lmu = log2(params$mu)
    modfit = loess(lsize ~ lmu,data=params)
    return(list(index=index,modfit=modfit,nbparams=params))
}

