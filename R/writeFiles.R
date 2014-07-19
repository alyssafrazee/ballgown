#' @title write files to disk from ballgown object
#'
#' @description create tablemaker-like files on disk from a ballgown object
#' @param gown ballgown object
#' @param dataDir top-level directory for sample-specific folders
#' @export
#' @examples \donttest{
#'   data(bg)
#'   writeFiles(bg, dataDir=getwd())
#'}
writeFiles = function(gown, dataDir){
    # write tablemaker output files to disk (from ballgown object)
    # I don't imagine this will be widely useful, but you never know.

    if(gown@meas != 'all'){
        msg = 'Can only write full ballgown objects to disk, for now.
            (You only have specific measurements in yours).'
        stop(.makepretty(msg))
    }

    sysoutdir = gsub(' ', '\\\\ ', dataDir) #vector, length 1
    for(f in sampleNames(gown)){
        sysout = paste0(sysoutdir, '/', f)
        if(.Platform$OS.type == 'windows'){
            shell(paste('mkdir ', sysout))
        }else{
            system(paste('mkdir -p', sysout))    
        }
        out = paste0(dataDir, '/', f)

        # i2t.ctab and e2t.ctab:
        write.table(indexes(gown)$i2t, file=paste0(out, '/i2t.ctab'), 
            quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')
        write.table(indexes(gown)$e2t, file=paste0(out, '/e2t.ctab'), 
            quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')

        # e_data.ctab
        cn = names(eexpr(gown, 'all'))
        e_data = eexpr(gown, 'all')[,c(1:5, which(grepl(f, cn)))]
        names(e_data)[6:ncol(e_data)] = c('rcount', 'ucount', 'mrcount',
            'cov', 'cov_sd', 'mcov', 'mcov_sd')
        write.table(e_data, file=paste0(out, '/e_data.ctab'), quote=FALSE, 
            col.names=TRUE, row.names=FALSE, sep='\t')

        # i_data.ctab
        cn = names(iexpr(gown, 'all'))
        i_data = iexpr(gown, 'all')[,c(1:5, which(grepl(f, cn)))]
        names(i_data)[6:ncol(i_data)] = c('rcount', 'ucount', 'mrcount')
        write.table(i_data, file=paste0(out, '/i_data.ctab'), quote=FALSE, 
            col.names=TRUE, row.names=FALSE, sep='\t')

        # t_data.ctab
        cn = names(texpr(gown, 'all'))
        t_data = texpr(gown, 'all')[,c(1:10, which(grepl(f, cn)))]
        names(t_data)[11:ncol(t_data)] = c('cov', 'FPKM')
        write.table(t_data, file=paste0(out, '/t_data.ctab'), quote=FALSE, 
            col.names=TRUE, row.names=FALSE, sep='\t')

    }
}
