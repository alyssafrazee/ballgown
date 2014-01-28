setMethod("texpr", "ballgown", function(x, meas="FPKM"){
  meas = match.arg(meas, c("cov","FPKM","all"))
  if(meas!="all"){
    expr = data(x)$trans[,-c(1:10)]
    expr = expr[,sapply(colnames(expr), function(x) strsplit(x,split="\\.")[[1]][1]==meas)]
    rownames(expr) = data(x)$trans$t_id
    expr = as.matrix(expr)
  }else{
    expr = data(x)$trans
  }
  return(expr)
})
