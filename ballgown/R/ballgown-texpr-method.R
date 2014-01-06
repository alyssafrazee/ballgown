setMethod("texpr", "ballgown", function(x, meas="all"){
  meas = match.arg(meas, c("cov","FPKM","all"))
  if(meas!="all"){
    expr = data(x)$trans[,-c(1:10)]
    expr = expr[,sapply(colnames(expr), function(x) strsplit(x,split="\\.")[[1]][1]==meas)]
    rownames(expr) = data(x)$trans$t_id
  }else{
    expr = data(x)$trans
  }
  return(expr)
})
