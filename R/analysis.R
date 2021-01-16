#' Variant calling using MuTECT2
#'
#' This function calls somatic variants in a pair of tumor-normal matched samples, or
#' just in a tumor sample if no matched sample is not available.
#'
#' @param tumor_bam Path to tumor bam file.
#' @param normal_bam Path to germline bam file.
#' @param bin_path Path to fastQC executable. Default path tools/gatk/gatk.
#' @param ref_genome Path to reference genome fasta file.
#' @param region Region to analyze. Optional
#' @param germ_resource Path to germline resources vcf file.
#' @param pon [Optional] Path to panel of normal.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export


vcf_mutect2=function(region="",bin_path="tools/gatk/gatk",tumor_bam="",normal_bam="",ref_genome="",germ_resource="",pon="",output_dir="",verbose=FALSE){
  sep="/"

  if(output_dir==""){
    sep=""
  }

  sample_name=ULPwgs::get_sample_name(tumor_bam[1])

  out_file=paste0(output_dir,sep,sample_name,"_MUTECT2_VARIANTS_VCF")
  if (!dir.exists(out_file)){
      dir.create(out_file)
  }
  reg=""
  if (region==""){
      out_file=paste0(out_file,"/",sample_name,".UNFILTERED_MUTECT2.vcf")
  }else{
      reg=paste0(" -L ",region)
      out_file=paste0(out_file,"/",sample_name,".",region,".UNFILTERED_MUTECT2.vcf")
  }



# TO DO FIX THIS MESS

  if (is.vector(tumor_bam)){
    tumor=paste0(" -I ",paste(tumor_bam,collapse=" -I "))
  }else{
    tumor=paste0(" -I ",tumor_bam)
  }
  norm=" "
  if (normal_bam!=""){
    if (is.vector(normal_bam)){
      norm=paste0(" -I ",paste(normal_bam,collapse=" -I ")," -normal ",paste(as.vector(sapply(normal_bam,FUN=ULPwgs::get_sample_name)),collapse=" -normal "))
  }else{
      norm=paste0(" -I ",normal_bam," -normal ",ULPwgs::get_sample_name(normal_bam))
      }
  }
  pn=" "
  if (pon!=""){
    pn=paste0(" --panel-of-normals ",pon)
  }
  if(verbose){
      print(paste0(bin_path," Mutect2 -R ",ref_genome,tumor, norm," --germline-resource ",germ_resource,pn," -O ",out_file,reg))

  }
  system(paste0(bin_path," Mutect2 -R ",ref_genome,tumor, norm, " --germline-resource ",germ_resource, pn, " -O ",out_file,reg))

}



#' VCF file concatenation
#'
#' This function concatenates VCF files found in a directory.
#'
#' @param bin_path Path to bcftools binary. Default tools/bcftools/bcftools.
#' @param vcf_dir Path to directory with vcf files to concatenate.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export



vcf_concatenate=function(bin_path="tools/bcftools/bcftools",vcf_dir="",verbose=FALSE,output_dir=""){
  sep="/"
  if(output_dir==""){
    sep=""
  }
  files=list.files(vcf_dir,full.names=TRUE)
  files=files[grepl(".vcf$",files)]
  sample_name=ULPwgs::get_sample_name(list.files(vcf_dir)[1])
  out_file_dir=paste0(output_dir,sep,sample_name,"_CONCATENATED")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  out_file=paste0(out_file_dir,"/",sample_name,".CONCATENATED.vcf")

  if(verbose){
    print(paste(bin_path,"concat -o",out_file, paste(files, collapse=' ' ) ))
  }
  system(paste(bin_path,"concat -o",out_file, paste(files, collapse=' ' )))
}

#' Predict variant effect using VEP
#'
#' This function predicts the effect of variant found in a VCF file
#'
#' @param bin_path Path to bcftools binary. Default tools/bcftools/bcftools.
#' @param vcf Path to vcf file.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export



vep=function(bin_path="tools/ensembl-vep/vep",vcf="",verbose=FALSE,output_dir=""){
  sep="/"
  if(output_dir==""){
    sep=""
  }

  sample_name=ULPwgs::get_sample_name(vcf)
  out_file_dir=paste0(output_dir,sep,sample_name,"_VEP")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  out_file=paste0(out_file_dir,"/",sample_name,".VEP.vcf")

  if(verbose){
    print(paste(bin_path,"-i",,vcf,"-o",out_file,"--cache --everything --force_overwrite --vcf"))
  }
  system(paste(bin_path,"-i",vcf,"-o",out_file,"--cache --everything --force_overwrite --vcf"))
}





#' VCF file concatenation
#'
#' This function merges VCF stat files found in a directory.
#'
#' @param bin_path Path to gatk binary. Default tools/gatk/gatk.
#' @param vcf_stats_dir Path to directory with vcf stats files to merge.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export



vcf_stats_merge=function(bin_path="tools/gatk/gatk",vcf_stats_dir="",verbose=FALSE,output_dir=""){
  sep="/"
  if(output_dir==""){
    sep=""
  }
  files=list.files(vcf_stats_dir,full.names=TRUE)
  files=files[grepl(".vcf.stats$",files)]
  sample_name=ULPwgs::get_sample_name(list.files(vcf_stats_dir)[1])
  out_file_dir=paste0(output_dir,sep,sample_name,"_MERGED_VCF_STATS")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  out_file=paste0(out_file_dir,"/",sample_name,".MERGED.vcf.stats")

  if(verbose){
    print(paste(bin_path,"MergeMutectStats -O",out_file," -stats ",paste(files, collapse=' -stats ' ) ))
  }
  system(paste(bin_path,"MergeMutectStats -O",out_file," -stats ", paste(files, collapse=' -stats ' )))
}





#' Sort VCF file
#'
#' This function takes a VCF file and sorts it genomic order
#'
#' @param bin_path [Required] Path to bcftools binary. Default tools/bcftools/bcftools.
#' @param vcf [Required] Path to VCF file to sort.
#' @param output_dir [Optional] Path to the output directory. Default working directory.
#' @param ram [Required] RAM in GB to allocate sorting. Default 4Gb.
#' @param verbose [Optional] Enables progress messages. Default False.
#' @param tmp_dir [Optional] Path to tmp file directory.
#' @export



vcf_sort=function(bin_path="tools/bcftools/bcftools",vcf="",verbose=FALSE,output_dir="",ram=4,tmp_dir=""){
  sep="/"
  if(output_dir==""){
    sep=""
  }

  sample_name=ULPwgs::get_sample_name(vcf)
  file_ext=ULPwgs::get_file_extension(vcf)
  out_file_dir=paste0(output_dir,sep,sample_name,"_SORTED.",toupper(file_ext))
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  out_file=paste0(out_file_dir,"/",sample_name)
  tmp=""
  if (!tmp_dir==""){
    tmp=paste0("-T ",tmp_dir)
  }

  if (verbose){
    print(paste0(bin_path," sort ",vcf," -m ",ram, "G ",tmp," -o ",out_file,".SORTED.",file_ext))
  }
  system(paste0(bin_path," sort ",vcf," -m ",ram, "G ",tmp," -o ",out_file,".SORTED.",file_ext))

}



#' Variant calling using MuTECT2 on parallel per genomic region
#'
#' This function calls somatic variants in a pair of tumor-normal matched samples, or
#' just in a tumor sample if no matched sample is not available.
#'
#' @param tumor_bam [Required] Path to tumor bam file.
#' @param normal_bam Path to germline bam file.
#' @param bin_path [Required] Path to gatk binary. Default path tools/gatk/gatk.
#' @param bin_path2 [Required] Path to bcftools binary. Default path tools/bcftools/bcftools.
#' @param bin_path3 [Required] Path to tabix binary. Default path tools/htslib/tabix.
#' @param bin_path4 [Required] Path to bgzip binary. Default path tools/htslib/bgzip.
#' @param ref_genome Path to reference genome fasta file.
#' @param germ_resource Path to germline resources vcf file.
#' @param pon [Optional] Path to panel of normal.
#' @param threads [Optional] Number of threads. Default 3
#' @param region_bed Path to bed file with regions to analyze.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export
#' @import pbapply


vcf_mutect2_parallel=function(bin_path="tools/gatk/gatk",bin_path2="tools/bcftools/bcftools",bin_path3="tools/htslib/bgzip",bin_path4="tools/htslib/tabix",tumor_bam="",normal_bam="",ref_genome="",germ_resource="",pon="",output_dir="",region_bed="",threads=3,verbose=FALSE){
  dat=read.table(region_bed)
  dat$V2=dat$V2+1
  dat=dat %>% dplyr::mutate(Region=paste0(sub("chr","",V1),":",V2,"-",V3))
  cl=parallel::makeCluster(threads)
  pbapply(X=dat[,c("Region"),drop=FALSE],1,FUN=vcf_mutect2,bin_path=bin_path,tumor_bam=tumor_bam,normal_bam=normal_bam,ref_genome=ref_genome,germ_resource=germ_resource,pon=pon,verbose=verbose,cl=cl)
  on.exit(parallel::stopCluster(cl))
  sep="/"
  if(output_dir==""){
    sep=""
  }
  sample_name=ULPwgs::get_sample_name(tumor_bam[1])
  out_file_dir=paste0(output_dir,sep,sample_name,"_MUTECT2_VARIANTS_VCF")
  vcf_concatenate(bin_path=bin_path2,vcf_dir=out_file_dir,output_dir=out_file_dir,verbose=verbose)
  vcf_sort(bin_path=bin_path2,vcf=paste0(out_file_dir,"/",sample_name,"_CONCATENATED","/",sample_name,".CONCATENATED.vcf.gz"),output_dir=out_file_dir,verbose=verbose)
  vcf_stats_merge(bin_path=bin_path,vcf_stats_dir=out_file_dir,output_dir=out_file_dir,verbose=verbose)
  vcf_filtering(bin_path=bin_path,bin_path2=bin_path3,bin_path3=bin_path4,ref_genome=ref_genome,unfil_vcf_stats=paste0(out_file_dir,"/",sample_name,"_SORTED.CONCATENATED.VCF","/",sample_name,".SORTED.CONCATENATED.vcf"),unfil_vcf_stats=paste0(out_file_dir,"/",sample_name,"_MERGED_VCF_STATS","/",sample_name,".MERGED.vcf.stats"),output_dir=out_file_dir,verbose=verbose)
}



#' VCF filtering using GATK
#'
#' This function filters VCF calls using GATK statistics
#'
#' @param bin_path Path to gatk binary. Default tools/gatk/gatk.
#' @param bin_path2 Path to bgzip binary. Default tools/htslib/bgzip.
#' @param bin_path3 Path to tabix binary. Default tools/htslib/tabix.
#' @param unfil_vcf Path to unfiltered vcf file.
#' @param unfil_vcf_stats Path to unfiltered vcf file stats.
#' @param ref_genome Path to reference genome fasta file.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export



vcf_filtering=function(bin_path="tools/gatk/gatk",bin_path2="tools/htslib/bgzip",bin_path3="tools/htslib/tabix",unfil_vcf="",ref_genome="",unfil_vcf_stats="",verbose=FALSE,output_dir=""){
  sep="/"
  if(output_dir==""){
    sep=""
  }
  sample_name=ULPwgs::get_sample_name(unfil_vcf)
  out_file_dir=paste0(output_dir,sep,sample_name,"_FILTERED")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  out_file=paste0(out_file_dir,"/",sample_name,".FILTERED.vcf")

  if(verbose){
    print(paste(bin_path,"FilterMutectCalls -O",out_file," -R ",ref_genome," -V ",unfil_vcf," -stats ",unfil_vcf_stats))
  }
  system(paste(bin_path,"FilterMutectCalls -O",out_file," -R ",ref_genome," -V ",unfil_vcf," -stats ",unfil_vcf_stats))
  bgzip(bin_path=bin_path2,file=out_file)
  tab_indx(bin_path=bin_path3,file=paste0(out_file,".gz"))
}




#' Variant calling using Platypus
#'
#' This function calls somatic variants in a pair of tumor-normal matched samples, or
#' just in a tumor sample if no matched sample is not available.
#'
#' @param tumor_bam Path to tumor bam file.
#' @param normal_bam Path to germline bam file.
#' @param bin_path Path to fastQC executable. Default path tools/platypus/Platypus.py.
#' @param ref_genome Path to reference genome fasta file.
#' @param vcf_overlay Path to vcf overlay to use as source.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export


vcf_platypus=function(bin_path="tools/platypus/Platypus.py",tumor_bam="",normal_bam="",ref_genome="",vcf_overlay="",output_dir="",verbose=FALSE){

  sep="/"

  if(output_dir==""){
    sep=""
  }

  sample_name=ULPwgs::get_sample_name(tumor_bam)

  out_file=paste0(output_dir,sep,sample_name,"_PLATYPUS_VARIANTS_VCF")
  if (!dir.exists(out_file)){
      dir.create(out_file)
  }
  if (is.vector(normal_bam)){
    norm=paste0(paste(normal_bam,collapse=","))
  }else{
    norm=normal_bam
  }

  if (is.vector(tumor_bam)){
    tumor=paste0(paste(tumor_bam,collapse=","))
  }else{
    tumor=tumor_bam
  }
  out_file=paste0(out_file,"/",sample_name,"_PLATYPUS.vcf")

# TO DO FIX THIS MESS

  if(verbose){
      print(paste0(bin_path," callVariants --refFile=",ref_genome, paste0(" --bamFiles=",tumor,",",norm)," --source=",vcf," --output=",out_file," --filterReadPairsWithSmallInserts=0 --minPosterior=0 --getVariantsFromBAMs=1 --logFileName=",paste0(out_file,".log")))

  }
  system(paste0(bin_path," callVariants --refFile=",ref_genome, paste0(" --bamFiles=",tumor,",",norm), " --source=",vcf," --output",out_file," --filterReadPairsWithSmallInserts=0 --minPosterior=0 --getVariantsFromBAMs=1 --logFileName=",paste0(out_file,".log")))





}
