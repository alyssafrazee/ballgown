setMethod("iexpr", "ballgown", function(x, meas="rcount"){
  meas = match.arg(meas, c("rcount","ucount","mrcount","all"))
  if(meas!="all"){
    expr = data(x)$intron[,-c(1:5)]
    expr = expr[,sapply(colnames(expr), function(x) strsplit(x,split="\\.")[[1]][1]==meas)]
    rownames(expr) = data(x)$intron$i_id
    expr = as.matrix(expr)
  }else{
    expr = data(x)$intron
  }
  return(expr)
})
