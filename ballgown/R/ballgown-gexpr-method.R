setMethod("gexpr", "ballgown", function(x){
    gnames = indexes(x)$t2g$g_id
    inds_by_gene = split(seq(along=gnames), gnames)
    tmeas = texpr(x, "FPKM")
    gid_by_exon = lapply(1:nrow(texpr(x)), function(i){rep(texpr(x)$gene_id[i], texpr(x)$num_exons[i])})
    ulstruct = unlist(structure(x)$trans)
    glist = split(ulstruct, unlist(gid_by_exon))
    glengths = sapply(width(reduce(glist)), sum)
    tlengths = sapply(width(structure(x)$trans), sum)
    tfrags = lapply(1:nrow(tmeas), function(i){
        (tlengths[i]/1000) * tmeas[i,]
    }) ## still a bit slow
    tfrags = matrix(unlist(tfrags, use.names=FALSE), nrow=length(tfrags), byrow=TRUE)
    expr = t(sapply(1:length(inds_by_gene), function(i){colSums(tfrags[inds_by_gene[[i]],,drop=FALSE]) / glengths[i]}))
    rownames(expr) = names(inds_by_gene)
    return(expr)
})

