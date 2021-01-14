#' Index tab separated genomic regions
#'
#' This function takes a .vcf.bgz tab separated genomic region file and generates an index for it
#'
#' @param bin_path [Required] Path to bcftools binary. Default tools/htslib/tabix.
#' @param file [Required] Path to VCF file to index.
#' @param verbose [Optional] Enables progress messages. Default False.
#' @export



tab_indx=function(bin_path="tools/htslib/tabix",file="",verbose=FALSE){

  if (verbose){
    print(paste0(bin_path," -fp vcf",file))
  }
  system(paste0(bin_path," -fp vcf ",file))
}



#' Bzips a VCF file
#'
#' This function takes a VCF file and Bgzips it
#'
#' @param bin_path [Required] Path to bcftools binary. Default tools/htslib/tabix.
#' @param file [Required] Path to VCF file to bgzip.
#' @param verbose [Optional] Enables progress messages. Default False.
#' @export



bgzip=function(bin_path="tools/htslib/bgzip",file="",verbose=FALSE){
  if (verbose){
    print(paste0(bin_path," -f ",file))
  }
  system(paste0(bin_path," -f ",file))
}
