#' extract exon-level expression measurements from ballgown objects
#' 
#' @name eexpr
#' @exportMethod eexpr
#' @docType methods
#' @rdname eexpr
#' @aliases eexpr,ballgown-method
#' @param x a ballgown object
#' @param meas type of measurement to extract. Can be "rcount", "ucount", "mrcount", "cov", "mcov",
#'             or "all". Default "rcount".
setMethod("eexpr", "ballgown", function(x, meas="rcount"){
    meas = match.arg(meas, c("rcount","ucount","mrcount","cov","mcov","all"))
    if(meas!="all"){
        expr = data(x)$exon[,-c(1:5)]
        expr = expr[,sapply(colnames(expr), function(x) strsplit(x,split="\\.")[[1]][1]==meas)]
        rownames(expr) = data(x)$exon$e_id
        expr = as.matrix(expr)
    }else{
        expr = data(x)$exon
    }
    return(expr)
})
