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

#' Splits multisample VCF file into multiple VCFs
#'
#' This function takes a multisample VCF and splits into multiple VCFS for each sample
#'
#' @param bin_path [Required] Path to bcftools binary. Default "tools/bcftools/bcftools"
#' @param file [Required] Path to VCF file to split.
#' @param verbose [Optional] Enables progress messages. Default False.
#' @param output_dir [Optional] Path to the output directory.
#' @export



split_vcf=function(bin_path="tools/bcftools/bcftools",vcf="",verbose=FALSE,output_dir=""){
    sep="/"

    if(output_dir==""){
      sep=""
    }

    sample_name=ULPwgs::get_sample_name(vcf)

    out_file=paste0(output_dir,sep,sample_name,"_SPLIT")
    if (!dir.exists(out_file)){
        dir.create(out_file)
    }

  if (verbose){
    print(paste(bin_path," query -l ",vcf))
    apply(names,1,function(x){print(paste(bin_path, "view -c1 -Oz -s",x, "-o", paste0(x,".vcf.gz"),vcf))
    })
  }

  names=system(paste(bin_path," query -l ",vcf), intern = TRUE, ignore.stderr = TRUE)
  apply(names,1,function(x){system(paste(bin_path, "view -c1 -Oz -s",x, "-o", paste0(out_file,"/",x,".vcf.gz"),vcf))
  })

}
