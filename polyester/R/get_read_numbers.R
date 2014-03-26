#' get read numbers from zero-inflated negative binomial simulation model
#' 
#' @param fc matrix of fold changes. Should have number of rows equal to the number of transcripts in the simulation, and number of columns equal to the number of samples in the simulation.
#' @param nbparams zero-inflated negative binomial parameters (estimated with \code{get_params})
#' @param modfit model defining the mean/dispersion relationship for transcript expression. Usually estimated with \code{get_params}.
#' @param libsizes how many reads should be generated for each sample? Should have length equal to number of columns of \code{fc}, i.e., \code{libsizes} is a vector with length equal to the number of samples in the simulation.
#' @param seed optional seed (for reproducibility)
#' @export
#' @author Jeff Leek

get_read_numbers = function(fc,nbparams,modfit,libsizes,seed=NULL){

    ## Set seed if asked for
    if(!is.null(seed)){set.seed(seed)}
  
    nsamples = dim(fc)[2]
    ntranscripts = dim(fc)[1]
  
    mumat = nbparams$mu * matrix(1,nrow=nrow(fc),ncol=ncol(fc))
    nreads = sizemat = mumat*NA
  
    for(i in 1:nrow(nreads)){
        sizemat[i,] = 2^predict(modfit,log2(mumat[i,]))
        nreads[i,] = rnbinom(dim(sizemat)[2],size=sizemat[i,],mu=mumat[i,])
        if(nbparams$p0[i] > 0){
            nreads[i,] = nreads[i,]*rbinom(nsamples,size=1,prob=(1-nbparams$p0[i]))
        }
    }
    nreads = nreads * fc
    nreads = round(scale(nreads,center=FALSE,scale=colSums(nreads)/libsizes)) 
    return(list(nreads=nreads,mumat=mumat,sizemat=sizemat))
}
