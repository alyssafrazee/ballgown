#' get feature-by-sample expression table from cuffdiff2/cummeRbund output
#' 
#' @param dbFile path to the database file from cummeRbund (available as an insilicodb download)
#' @param feature one of \code{"isoform"}, \code{"gene"}, \code{"CDS"}, or \code{"TSS"} -- what genomic feature do you want expression measurements for? Default \code{"isoform"}
#' @param meas one of \code{"fpkm"}, \code{"raw_frags"}, \code{"internal_scaled_frags"}, \code{"external_scaled_frags"} -- which measurement from the Cuffdiff output table do you want? Default \code{"fpkm"}.
#' @return matrix with \code{feature}s in rows, samples (replicates) in columns, and \code{meas}urements in cells
#' @details on a test dataset with 37K transcripts and 24 reps, this took about 30 seconds
#' @author Alyssa Frazee
#' @export

getInsdbTable = function(dbFile, feature="isoform", meas="fpkm"){
    require(cummeRbund)
    require(reshape)

    feature = match.arg(feature, c("isoform", "gene", "CDS", "TSS"))
    meas = match.arg(meas, c("fpkm", "raw_frags", "internal_scaled_frags", "external_scaled_frags"))

    cuff = readCufflinks(dbFile=dbFile)

    if(feature == "isoform"){
        data_tmp = cast(repFpkm(isoforms(cuff)), isoform_id ~ rep_name, value=meas)
    }else if(feature == "gene"){
        data_tmp = cast(repFpkm(genes(cuff)), gene_id ~ rep_name, value=meas)
    }else if(feature == "CDS"){
        data_tmp = cast(repFpkm(CDS(cuff)), CDS_id ~ rep_name, value=meas)
    }else{
        data_tmp = cast(repFpkm(TSS(cuff)), TSS_id ~ rep_name, value=meas)
    }

    final_data = as.matrix(data_tmp[,-1])
    rownames(final_data) = data_tmp[,1]
    colnames(final_data) = names(data_tmp)[-1]

    return(final_data)

}

isoform_fpkm = getInsdbTable("GSE37764GPL10999_DGE_a9dc2c94672e4a51c036c76be9508164.db")

