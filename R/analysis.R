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
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @param pon [Optional] Path to panel of normal.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export


vcf_mutect2=function(region="",bin_path="tools/gatk/gatk",tumor_bam="",normal_bam="",ref_genome="",germ_resource="",pon="",output_dir="",output_name="",verbose=FALSE){
  sep="/"

  if(output_dir==""){
    sep=""
  }

  if (output_name==""){
    sample_name=ULPwgs::get_sample_name(tumor_bam[1])
  }else{
    sample_name=output_name
  }

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

#' Variant calling using bcftools
#'
#' This function calls variants in a sample.
#'
#'
#' @param bam Path to bam file.
#' @param bin_path Path to fastQC executable. Default path tools/gatk/gatk.
#' @param ref_genome Path to reference genome fasta file.
#' @param region Region to analyze. Optional
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.


vcf_bcftools=function(region="",bin_path="tools/bcftools/bcftools",bam="",ref_genome="",output_dir="",verbose=FALSE){

  sep="/"

  if(output_dir==""){
    sep=""
  }

  sample_name=ULPwgs::get_sample_name(bam[1])

  out_file=paste0(output_dir,sep,sample_name,"_BCF_MPILEUP_VARIANTS_VCF")
  if (!dir.exists(out_file)){
      dir.create(out_file)
  }
  reg=""
  if (region==""){
      out_file=paste0(out_file,"/",sample_name,".UNFILTERED_MPILEUP.vcf")
  }else{
      reg=paste0(" -r ",region)
      out_file=paste0(out_file,"/",sample_name,".",region,".UNFILTERED_MPILEUP.vcf")
  }

  if (is.vector(bam)){
    bam=paste(bam,collapse=" ")
  }else{
    bam=paste0(" ",bam)
  }

  if(verbose){
      print(paste0(bin_path," mpileup ",reg," -f ",paste(ref_genome, bam)," | " ,bin_path," call --multiallelic-caller --variants-only -Ov >",out_file))

  }
  system(paste0(bin_path," mpileup ", reg ," -f ",paste(ref_genome, bam)," | " ,bin_path," call --multiallelic-caller --variants-only -Ov >",out_file))

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
#' @param threads Number of threads to use. Default 3.
#' @export



call_vep=function(bin_path="tools/ensembl-vep/vep",vcf="",verbose=FALSE,output_dir="",threads=3){
  sep="/"
  if(output_dir==""){
    sep=""
  }

  sample_name=ULPwgs::get_sample_name(vcf)
  file_ext=ULPwgs::get_file_extension(vcf)
  out_file_dir=paste0(output_dir,sep,sample_name,"_VEP")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  out_file=paste0(out_file_dir,"/",sample_name,".VEP.",file_ext)

  if(verbose){
    print(paste(bin_path,"-i",,vcf,"-o",out_file,"--cache --port 3337 --everything --force_overwrite --vcf --fork",threads))
  }
  system(paste(bin_path,"-i",vcf,"-o",out_file,"--cache --port 3337 --everything --force_overwrite --vcf --fork",threads))
}



call_clonet=function(bin_path="tools/CLONET/CLONET",vcf="",verbose=FALSE,output_dir="",threads=3){
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
    print(paste(bin_path,"-i",,vcf,"-o",out_file,"--cache --port 3337 --everything --force_overwrite --vcf --fork",threads))
  }
  system(paste(bin_path,"-i",vcf,"-o",out_file,"--cache --port 3337 --everything --force_overwrite --vcf --fork",threads))
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
  sample_name=ULPwgs::get_sample_name(files[1])
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
  out_file_dir=paste0(output_dir,sep,sample_name,"_SORTED")
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
#' @param tumor_bam [REQUIRED] Path to tumor bam file/s.
#' @param normal_bam [OPTIONAL] Path to germline bam file/s.
#' @param bin_path [REQUIRED] Path to gatk binary. Default path tools/gatk/gatk.
#' @param bin_path2 [REQUIRED] Path to bcftools binary. Default path tools/bcftools/bcftools.
#' @param bin_path3 [REQUIRED] Path to tabix binary. Default path tools/htslib/tabix.
#' @param bin_path4 [REQUIRED] Path to bgzip binary. Default path tools/htslib/bgzip.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param germ_resource [REQUIRED] Path to germline resources vcf file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @param pon [OPTIONAL] Path to panel of normal.
#' @param threads [OPTIONAL] Number of threads. Default 3
#' @param chr_filter [OPTIONAL] Chromosomes to analyze. canonical/autosomal/all or a list of chromosomes
#' @param region_bed [REQUIRED] Path to bed file with regions to analyze.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @export
#' @import pbapply


call_mutect2_parallel=function(bin_path="tools/gatk/gatk",bin_path2="tools/bcftools/bcftools",bin_path3="tools/htslib/bgzip",bin_path4="tools/htslib/tabix",tumor_bam="",normal_bam="",ref_genome="",germ_resource="",pon="",output_dir="",output_name="",region_bed="",threads=3,verbose=FALSE,chr_filter="canonical"){

  if (output_name==""){
    sample_name=ULPwgs::get_sample_name(tumor_bam[1])
  }else{
    sample_name=output_name
  }
  dat=read.table(region_bed)

  ## Dont include unplaced contigs, a.k.a only autosomal+sexual+MT chromosomes
  if(chr_filter=="canonical"){
    dat=dat[grepl(paste0(paste0(paste0("^",c(1:22,"X","Y","M","MT")),"$"),collapse="|"),dat$V1),]
  ## Include autosomal chromosomal and nothing else
  }else if(chr_filter=="autosomal"){
    dat=dat[grepl(paste0(paste0(paste0("^",c(1:22)),"$"),collapse="|"),dat$V1),]
  ## Include only specific chromosomes
  }else if(chr_filter!="all"){
    dat=dat[grepl(paste0(paste0(paste0("^",chr_filter),"$"),collapse="|"),dat$V1),]
  }
  dat$V2=dat$V2+1
  dat=dat %>% dplyr::mutate(Region=paste0(sub("chr","",V1),":",V2,"-",V3))
  cl=parallel::makeCluster(round(threads/4), digits = 0)
  pbapply(X=dat[,c("Region"),drop=FALSE],1,FUN=vcf_mutect2,bin_path=bin_path,tumor_bam=tumor_bam,normal_bam=normal_bam,ref_genome=ref_genome,germ_resource=germ_resource,pon=pon,output_dir=output_dir,output_name=output_name,verbose=verbose,cl=cl)
  on.exit(parallel::stopCluster(cl))
  sep="/"
  if(output_dir==""){
    sep=""
  }
  out_file_dir=paste0(output_dir,sep,sample_name,"_MUTECT2_VARIANTS_VCF")
  vcf_concatenate(bin_path=bin_path2,vcf_dir=out_file_dir,output_dir=out_file_dir,verbose=verbose)
  vcf_sort(bin_path=bin_path2,vcf=paste0(out_file_dir,"/",sample_name,"_CONCATENATED","/",sample_name,".CONCATENATED.vcf.gz"),output_dir=out_file_dir,verbose=verbose)
  vcf_stats_merge(bin_path=bin_path,vcf_stats_dir=out_file_dir,output_dir=out_file_dir,verbose=verbose)
  vcf_filtering(bin_path=bin_path,bin_path2=bin_path3,bin_path3=bin_path4,ref_genome=ref_genome,unfil_vcf=paste0(out_file_dir,"/",sample_name,"_SORTED","/",sample_name,".SORTED.CONCATENATED.vcf.gz"),
  unfil_vcf_stats=paste0(out_file_dir,"/",sample_name,"_MERGED_VCF_STATS","/",sample_name,".MERGED.vcf.stats"),output_dir=out_file_dir,verbose=verbose)
  ##system(paste0("rm -rf ",out_file_dir,"/*:*"))
  ##system(paste0("rm -rf ",out_file_dir,"/",sample_name,"_CONCATENATED"))
  ##system(paste0("rm -rf ",out_file_dir,"/",sample_name,"_MERGED_VCF_STATS"))
  ##system(paste0("rm -rf ",out_file_dir,"/",sample_name,"_SORTED.CONCATENATED.VCF.GZ"))


}

#' Variant calling using bcftools on parallel per genomic region
#'
#' This function calls SNPs in a sample
#'
#' @param bam [Required] Path to bam file.
#' @param bin_path [Required] Path to bcftools binary. Default path tools/bcftools/bcftools.
#' @param ref_genome Path to reference genome fasta file.
#' @param threads [Optional] Number of threads. Default 3
#' @param region_bed Path to bed file with regions to analyze.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export
#' @import pbapply

call_bcftools_parallel=function(bin_path="tools/bcftools/bcftools",bam="",ref_genome="",output_dir="",region_bed="",threads=3,verbose=FALSE){
  dat=read.table(region_bed)
  dat$V2=dat$V2+1
  dat=dat %>% dplyr::mutate(Region=paste0(sub("chr","",V1),":",V2,"-",V3))
  dat=dat %>% dplyr::filter(!grepl("_|MT|M",Region))
  cl=parallel::makeCluster(threads)
  pbapply(X=dat[,c("Region"),drop=FALSE],1,FUN=vcf_bcftools,bin_path=bin_path,bam=bam,ref_genome=ref_genome,verbose=verbose,cl=cl)
  on.exit(parallel::stopCluster(cl))
  sep="/"
  if(output_dir==""){
    sep=""
  }
  sample_name=ULPwgs::get_sample_name(bam[1])
  out_file_dir=paste0(output_dir,sep,sample_name,"_BCF_MPILEUP_VARIANTS_VCF")
  vcf_concatenate(bin_path=bin_path,vcf_dir=out_file_dir,output_dir=out_file_dir,verbose=verbose)
  vcf_sort(bin_path=bin_path,vcf=paste0(out_file_dir,"/",sample_name,"_CONCATENATED","/",sample_name,".CONCATENATED.vcf"),output_dir=out_file_dir,verbose=verbose)
  system(paste("cp", out_file, paste0(out_file,".tmp")))
  bgzip(bin_path=bin_path2,file=out_file)
  tab_indx(bin_path=bin_path3,file=paste0(out_file,".gz"))
  system(paste("cp", paste0(out_file,".tmp"), out_file))
  system(paste("rm -rf", paste0(out_file,".tmp")))
}

#' Variant calling using bcftools on parallel per genomic region
#'
#' This function calls SNPs in a sample
#'
#' @param bam [Required] Path to bam file.
#' @param bin_path [Required] Path to bcftools binary. Default path tools/bcftools/bcftools.
#' @param ref_genome Path to reference genome fasta file.
#' @param threads [Optional] Number of threads. Default 3
#' @param region_bed Path to bed file with regions to analyze.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export
#' @import pbapply

call_bcftools_parallel=function(bin_path="tools/bcftools/bcftools",bam="",ref_genome="",output_dir="",region_bed="",threads=3,verbose=FALSE){
  dat=read.table(region_bed)
  dat$V2=dat$V2+1
  dat=dat %>% dplyr::mutate(Region=paste0(sub("chr","",V1),":",V2,"-",V3))
  dat=dat %>% dplyr::filter(!grepl("_|MT|M",Region))
  cl=parallel::makeCluster(threads)
  pbapply(X=dat[,c("Region"),drop=FALSE],1,FUN=vcf_bcftools,bin_path=bin_path,bam=bam,ref_genome=ref_genome,verbose=verbose,cl=cl)
  on.exit(parallel::stopCluster(cl))
  sep="/"
  if(output_dir==""){
    sep=""
  }
  sample_name=ULPwgs::get_sample_name(bam[1])
  out_file_dir=paste0(output_dir,sep,sample_name,"_BCF_MPILEUP_VARIANTS_VCF")
  vcf_concatenate(bin_path=bin_path,vcf_dir=out_file_dir,output_dir=out_file_dir,verbose=verbose)
  vcf_sort(bin_path=bin_path,vcf=paste0(out_file_dir,"/",sample_name,"_CONCATENATED","/",sample_name,".CONCATENATED.vcf"),output_dir=out_file_dir,verbose=verbose)
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
  system(paste("cp", out_file, paste0(out_file,".tmp")))
  bgzip(bin_path=bin_path2,file=out_file)
  tab_indx(bin_path=bin_path3,file=paste0(out_file,".gz"))
  system(paste("cp", paste0(out_file,".tmp"), out_file))
  system(paste("rm -rf", paste0(out_file,".tmp")))
}



#' VCF filtering using bcftools
#'
#' This function filters VCF calls using bcftools
#'
#' @param bin_path Path to gatk binary. Default tools/gatk/gatk.
#' @param bin_path2 Path to bgzip binary. Default tools/htslib/bgzip.
#' @param bin_path3 Path to tabix binary. Default tools/htslib/tabix.
#' @param unfil_vcf Path to unfiltered vcf file.
#' @param qual Quality filter. Default 30.
#' @param mq Mapping quality filter. Default 40.
#' @param ref [Optional] Filter variants with ID.
#' @param state [Optional] Variant state to select. Options: het/homo
#' @param type [Optional] Variant type to include. Options: snp/indel. Default includes both
#' @param filter [Optional] Filter to include Options: PASS/.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export



vcf_filter_variants=function(unfil_vcf="",bin_path="tools/bcftools/bcftools",bin_path2="tools/htslib/bgzip",bin_path3="tools/htslib/tabix",qual=30,mq=40,state="",ref="",type="",filter="",verbose=FALSE,output_dir=""){

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


  if (state!=""){
    state=paste0(" & GT[0]=\"",state,"\" ")
  }

  if (ref!=""){
    ref=paste0(" & ID!=\"",ref,"\" ")
  }

  if (type!=""){
    type=paste0(" & TYPE=\"",type,"\" ")
  }

  if (filter!=""){
    filter=paste0(" & FILTER=\"",filter,"\" ")
  }


  if(verbose){
    print(paste(bin_path,"view  -i \'%QUAL>",qual,state,type,,filter,ref,"& MQ>",mq,"\'",unfil_vcf,">",out_file))
  }
  system(paste(bin_path,"view  -i \'%QUAL>",qual,state,type,filter,ref,"& MQ>",mq,"\'",unfil_vcf,">",out_file))
  system(paste("cp", out_file, paste0(out_file,".tmp")))
  bgzip(bin_path=bin_path2,file=out_file)
  tab_indx(bin_path=bin_path3,file=paste0(out_file,".gz"))
  system(paste("cp", paste0(out_file,".tmp"), out_file))
  system(paste("rm -rf", paste0(out_file,".tmp")))
}

#' This function filter and formats heterozygous SNP data for downstream analysis using clonet
#'
#' This function takes a path of a directory of unfiltered VCFs files generated using Platypus and
#' the a path of the directory of BAMs from which those VCFs were generated.
#'
#' @param bin_path Path to gatk binary. Default tools/bcftools/bcftools.
#' @param bin_path2 Path to bgzip binary. Default tools/htslib/bgzip.
#' @param bin_path3 Path to tabix binary. Default tools/htslib/tabix.
#' @param bin_path4 Path to ASEQ binary. Default tools/ASEQ/binaries/linux64/ASEQ
#' @param unfil_vcf_dir Path to unfiltered VCF file directory.
#' @param qual Quality filter. Default 30.
#' @param mq Mapping quality filter. Default 40.
#' @param bam_dir Path to BAM file directory.
#' @param patient_id Patient identifier.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads to use. Default 3
#' @export

format_SNP_data=function(bin_path="tools/bcftools/bcftools",bin_path2="tools/htslib/bgzip",bin_path3="tools/htslib/tabix",bin_path4="tools/ASEQ/binaries/linux64/ASEQ",unfil_vcf_dir="",bam_dir="",qual=30,mq=40,patient_id="",verbose=FALSE,output_dir="",threads=3){
  sep="/"
  if(output_dir==""){
    sep=""
  }
  out_file_dir=paste0(output_dir,sep,patient_id,"_PROCESSED_SNPs")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }
  files0=list.files(unfil_vcf_dir,recursive=TRUE,full.names=TRUE,pattern=patient_id)
  files0=files0[grepl("vcf.gz$",files0)]

  cl=parallel::makeCluster(threads)
  pbapply::pbapply(X=as.data.frame(files0),1,FUN=vcf_filter_variants,bin_path=bin_path,bin_path2=bin_path2,bin_path3=bin_path3,qual=qual,mq=mq,state="het",type="snp",filter="PASS",verbose=verbose,output_dir=out_file_dir,cl=cl)
  files1=list.files(out_file_dir,recursive=TRUE,full.names=TRUE)
  files1=files1[grepl("vcf$",files1)]

  pbapply::pbapply(X=as.data.frame(files1),1,FUN=vcf_format,bin_path=bin_path,bin_path2=bin_path2,bin_path3=bin_path3,expr="'%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO\\n'",verbose=verbose,output_dir=out_file_dir,cl=cl)
  files2=list.files(out_file_dir,recursive=TRUE,full.names=TRUE,pattern="FORMATED")
  files2=files2[grepl("vcf$",files2)]
  files2=as.data.frame(files2)
  names(files2)="VCF_path"
  files2$Sample=lapply(files2$VCF_path,FUN=ULPwgs::get_sample_name)

  files3=list.files(bam_dir,recursive=TRUE,full.names=TRUE,pattern=patient_id)
  files3=files3[grepl("bam$",file3)]
  files3=as.data.frame(files3)
  names(files3)="BAM_path"
  files3$Sample=lapply(files3$BAM_path,FUN=ULPwgs::get_sample_name)
  files=left_join(files2,files3,by="Sample")
  pbapply::pbapply(FUN=function(x){
    call_ASEQ(vcf=files[x,]$VCF_path,bin_path=bin_path,bam=files[x,]$BAM_path,output_dir=output_dir,threads=1,verbose=verbose)},cl=cl)

  files3=list.files(bam_dir,recursive=TRUE,full.names=TRUE,pattern="PILEUP.ASEQ")
  pbapply::pbapply(X=as.data.frame(files3),1,FUN=format_ASEQ_pileup,verbose=verbose,output_dir=out_file_dir,cl=cl)
  on.exit(parallel::stopCluster(cl))
}

#' VCF annotation using bcftools
#'
#' This function annotates VCF file using bcftools
#'
#' @param bin_path Path to gatk binary. Default tools/gatk/gatk.
#' @param bin_path2 Path to bgzip binary. Default tools/htslib/bgzip.
#' @param bin_path3 Path to tabix binary. Default tools/htslib/tabix.
#' @param vcf Path to vcf file.
#' @param db Path to vcf file with SNPs database.
#' @param state Variant state. Default het.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Additional in/out compression threads. Default 1.
#' @export



vcf_annotate=function(bin_path="tools/bcftools/bcftools",bin_path2="tools/htslib/bgzip",bin_path3="tools/htslib/tabix",vcf="",db="",verbose=FALSE,output_dir="",threads=1){
  sep="/"
  if(output_dir==""){
    sep=""
  }
  sample_name=ULPwgs::get_sample_name(vcf)
  out_file_dir=paste0(output_dir,sep,sample_name,"_ANNOTATED")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  out_file=paste0(out_file_dir,"/",sample_name,".ANNOTATED.vcf")

  if(verbose){
    print(paste(bin_path,"annotate -a ",db,"-c ID --threads",threads, vcf,">",out_file))
  }
  system(paste(bin_path,"annotate -a ",db,"-c ID --threads", threads,vcf,">",out_file))
  system(paste("cp", out_file, paste0(out_file,".tmp")))
  bgzip(bin_path=bin_path2,file=out_file)
  tab_indx(bin_path=bin_path3,file=paste0(out_file,".gz"))
  system(paste("cp", paste0(out_file,".tmp"), out_file))
  system(paste("rm -rf", paste0(out_file,".tmp")))
}


#' Variant calling of Somatic and Germline (MuTECT2 and Platypus) mutations and their subsequent annotation (VEP)
#'
#' This function calls somatic and germline mutations, in previously pre-processed sequencing data.
#' This function takes a path to a directory with BAM files and a patients ID. Then it subsets all BAM
#' files specific to the patient, identifying between cancer and normal Samples using the germline identifier.
#' Afterwards, it calls the MuTECT2 somatic variant caller, to indentify somatics variants, and Platypus genotyper all variants found in the samples.
#' Finally, it calls Varant Effect Predictor (VEP) and annotates all the variants.
#'
#'
#' @param bin_path [REQUIRED] Path to gatk binary. Default tools/gatk/gatk.
#' @param bin_path2 [REQUIRED] Path to bcftools binary. Default tools/bcftools/bcftools.
#' @param bin_path3 [REQUIRED] Path to bgzip binary. Default tools/htslib/bgzip.
#' @param bin_path4 [REQUIRED] Path to tabix binary. Default tools/htslib/tabix.
#' @param bin_path5 [REQUIRED] Path to platypus binary. Default tools/platypus/Platypus.py.
#' @param bin_path6 [REQUIRED] Path to vep binary. tools/ensembl-vep/vep
#' @param bam_dir [REQUIRED] Path to directory with BAM files.
#' @param patient_id [REQUIRED] Patient ID to analyze. Has to be in file names to subselect samples.
#' @param germ_pattern [REQUIRED] Pattern used to identify germline samples. Ex GL
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param germ_resource [REQUIRED] Path to germline resources vcf file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @param pon [OPTIONAL] Path to panel of normal.
#' @param threads [OPTIONAL] Number of threads. Default 3
#' @param region_bed [REQUIRED] Path to bed file with regions to analyze.
#' @param chr_filter [OPTIONAL] Chromosomes to analyze. canonical/autosomal/all or a list of chromosomes
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @export

call_variants=function(bin_path="tools/gatk/gatk",bin_path2="tools/bcftools/bcftools",bin_path3="tools/htslib/bgzip",bin_path4="tools/htslib/tabix",bin_path5="tools/platypus/Platypus.py",bin_path6="tools/ensembl-vep/vep",bam_dir="",patient_id="",germ_pattern="GL",ref_genome="",germ_resource="",pon="",output_dir="",region_bed="",chr_filter="canonical",threads=3,verbose=FALSE){
    sep="/"
    if(output_dir==""){
      sep=""
    }
    out_file_dir=paste0(output_dir,sep,patient_id,"_VARIANTS")
    if (!dir.exists(out_file_dir)){
        dir.create(out_file_dir)
    }
    files=list.files(bam_dir,recursive=TRUE,full.names=TRUE,pattern=patient_id)
    files=files[grepl("bam$",files)]
    tumor_bam=files[!grepl(germ_pattern,files)]
    normal_bam=files[grepl(germ_pattern,files)]
    call_mutect2_parallel(bin_path=bin_path,bin_path2=bin_path2,bin_path3=bin_path3,bin_path4=bin_path4,tumor_bam=tumor_bam,normal_bam=normal_bam,ref_genome=ref_genome,germ_resource=germ_resource,pon=pon,output_dir=out_file_dir,region_bed=region_bed,threads=threads,verbose=verbose,output_name=patient_id,chr_filter=chr_filter)
    call_platypus(bin_path=bin_path5,bin_path2=bin_path3,bin_path3=bin_path4,tumor_bam=tumor_bam,normal_bam=normal_bam,ref_genome=ref_genome,vcf_overlay=paste0(out_file_dir,"/",patient_id,"_MUTECT2_VARIANTS_VCF/",patient_id,"_FILTERED/",patient_id,".FILTERED.vcf"),output_dir=out_file_dir,verbose=verbose,threads=threads,output_name=patient_id)
    call_vep(bin_path=bin_path6,vcf=paste0(out_file_dir,"/",patient_id,"_MUTECT2_VARIANTS_VCF/",patient_id,"_FILTERED/",patient_id,".FILTERED.vcf"),verbose=verbose,output_dir=paste0(out_file_dir,"/",patient_id,"_MUTECT2_VARIANTS_VCF"),threads=threads)
    call_vep(bin_path=bin_path6,vcf=paste0(out_file_dir,"/",patient_id,"_PLATYPUS_VARIANTS_VCF/",patient_id,".PLATYPUS.vcf"),verbose=verbose,output_dir=paste0(out_file_dir,"/",patient_id,"_PLATYPUS_VARIANTS_VCF"),threads=threads)
}


#' VCF formating using bcftools
#'
#' This function formats VCF file using bcftools
#'
#' @param bin_path Path to gatk binary. Default tools/gatk/gatk.
#' @param bin_path2 Path to bgzip binary. Default tools/htslib/bgzip.
#' @param bin_path3 Path to tabix binary. Default tools/htslib/tabix.
#' @param vcf Path to vcf file.
#' @param expr Expresion by which to format VCF. Default %CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO\\n
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export

vcf_format=function(vcf="",bin_path="tools/bcftools/bcftools",bin_path2="tools/htslib/bgzip",bin_path3="tools/htslib/tabix",expr="'%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO\\n'",verbose=FALSE,output_dir=""){
  sep="/"
  if(output_dir==""){
    sep=""
  }
  sample_name=ULPwgs::get_sample_name(vcf)
  out_file_dir=paste0(output_dir,sep,sample_name,"_FORMATED")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  out_file=paste0(out_file_dir,"/",sample_name,".FORMATED.vcf")

  if(verbose){
    print(paste(bin_path,"query -f ",expr, vcf,">",out_file))
  }
  system(paste(bin_path,"query -f ",expr, vcf,">",out_file))
  system(paste("cp", out_file, paste0(out_file,".tmp")))
  bgzip(bin_path=bin_path2,file=out_file)
  tab_indx(bin_path=bin_path3,file=paste0(out_file,".gz"))
  system(paste("cp", paste0(out_file,".tmp"), out_file))
  system(paste("rm -rf", paste0(out_file,".tmp")))
}


#' ASEQ pileup formatting for downstream analysis

#' This function formats ASEQ pileup file for downstream analysis.
#'
#' @param file Path to file to format.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @import data.table
#' @export


format_ASEQ_pileup=function(file="",verbose=FALSE,output_dir=""){
  sep="/"
  if(output_dir==""){
    sep=""
  }
  sample_name=ULPwgs::get_sample_name(file)
  out_file_dir=paste0(output_dir,sep,sample_name,"_FORMATTED_ASEQ_PILEUP")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  out_file=paste0(out_file_dir,"/",sample_name,".snps")
  dat=read.table(file,header=TRUE)
  dat=dat %>% dplyr::select(chr,pos,dbsnp,ref,alt,A,C,G,T,RD)
  dat=dat %>% dplyr::filter(!nchar(as.character(alt))>1)
  dat=as.data.table(dat)
  dat=dat[,Value := get(as.character(alt)), by = alt]
  dat=dat %>% dplyr::mutate(af=Value/RD,cov=RD) %>% dplyr::select(chr, pos,dbsnp,ref,alt,A,C,G,T,af,cov) %>% dplyr::rename (rsid="dbsnp")
  write.table(dat,file=out_file,quote=FALSE,row.names=FALSE,sep="\t")
}





#' Call segments for panel, exome and WGS data using CNVkit

#' This function takes a single or multiple tumor samples and a single or multiple normal samples,
#' creates a reference of a pool of normal samples, estimates the coverage and normalizes its log2 ratio using gc and accessibility information.
#' Ultimately, it calls the segments of binned data and identifies the state of CNA.
#'
#' @param bin_path Path to cnvkit binary.
#' @param tumor_samples [REQUIRED] Path to tumor samples bam files. For multiple samples pass as vector with their paths.
#' @param normal_samples [OPTIONAL] Path to normal samples bam files. For multiple samples pass as vector with their paths. Used to create a reference.
#' @param targets [OPTIONAL] Path to exome capture targets. Required if targeted or panel.
#' @param access [OPTIONAL] Path to bed with accessibility information.
#' @param diagram [DEFAULT==TRUE] Plot diagram of segments.
#' @param scatter [DEFAULT==TRUE] Plot scatter plot of segments.
#' @param ref_output [OPTIONAL] Name of the reference file output.
#' @param threads [DEFAULT==3] Number of threads to use.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param male [DEFAULT==TRUE] Reference is a male.
#' @param verbose [DEFAULT==FALSE] Enables progress messages.
#' @export


call_segments=function(bin_path="~/tools/cnvkit/cnvkit.py",tumor_samples="",normal_samples="",targets="",fasta="",access="",ref_output="",output_dir="",diagram=TRUE,scatter=TRUE,threads=3,verbose=FALSE,male=TRUE){


  if (is.vector(tumor_samples) & length(tumor_samples)>1){
    tumor_samples=paste(tumor_samples,collapse=" ")
  }else{
    tumor_samples=paste0(" ",tumor_samples)
  }


  if (is.vector(tumor_samples) & length(tumor_samples)>1){
    normal_samples=paste(normal_samples,collapse=" ")
  }else{
    normal_samples=paste0(" --normal ",normal_samples)
  }

  if(male){
    mal=" -y "
  }
  add=""
  if(scatter){
    add=paste(add," --scatter ")
  }
  if(diagram){
    add=paste(add," --diagram ")
  }
  if (output_dir!=""){
    output_dir=paste(" --output-dir ",output_dir)
  }

  if (ref_output!=""){
    ref_output=paste(" --output-reference ",ref_output)
  }

  if (access!=""){
    access=paste(" --access",access)
  }


  if (fasta!=""){
    fasta=paste(" --fasta",fasta)
  }

  if (targets!=""){
    targets=paste(  " --targets ",targets)
  }

  if(verbose){
    print(paste(bin_path,"batch ",tumor_samples,normal_samples,targets,fasta,mal,ref_output,output_dir," --p ",threads,add))
  }
  system(paste(bin_path,"batch ",tumor_samples,normal_samples,targets,fasta,mal,ref_output,output_dir," --p ",threads,add))

}

#' Format segmentation data for downstream analysis using CLONETv2
#'
#' This function takes the segmentation data of multiple samples produced by any segment caller (mainly CNVkit) and generates
#' a single BED file with all the segment information (chr/start/end/log2) with an additional column for sample ID.
#' The default columns to select from original segmentation files are 1,2,3,5, which correspond to chr/start/end/log2 in CNVkit.
#' Use argument cols_to_keep to select other columns if needed or if order is different from different segment callers.
#'
#' @param seg_file [REQUIRED] Path/s to segmentation file/s. dir_segment and seg_file are mutually excluding.
#' @param dir_segment [REQUIRED] Path to directory with segmentation files.
#' @param pattern [FALSE] Pattern to use if directory for segmentation files is given.
#' @param cols_to_keep [DEFAULT==c(1,2,3,5)] Columns to keep from original segmentation bed files
#' @param output_dir Path to the output directory.
#' @param output_name Name of the file to output.
#' @param verbose Enables progress messages. Default False.
#' @export


format_segment_data=function(seg_file="",dir_segment="",pattern="",cols_to_keep=c(1,2,3,5),output_dir="",output_name="",verbose=FALSE){

  if(dir_segment!="" & seg_file!=""){

    stop("Only seg_file or dir_segment can be provided, not both.")
  }

  if (dir_segment!=""){
    files=list.files(path=dir_segment,pattern=pattern,full.names=TRUE)
  }else {
    files=seg_file
  }



  data=lapply(files,FUN=function(x) {dat=read.table(file=x,header=TRUE);
  dat$sample=ULPwgs::get_sample_name(x);
  return (dat)
  })
  data=dplyr::bind_rows(data)
  data=data[,c(cols_to_keep,ncol(data))]
  names(data)=c("chr","start","end","log2","sample")
  sep="/"
  if(output_dir==""){
    sep=""
  }

  if(output_name==""){
    output_name="SegmentData"
  }
  out_file_dir=paste0(output_dir,sep,"FORMATTED_SEGMENTS")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  out_file=paste0(out_file_dir,"/",output_name,".bed")
  write.table(data,file=out_file,quote=FALSE,row.names=FALSE,sep="\t")
}



#' Wrapper around ASEQ tool for pileup data.
#'
#'
#' This function takes a VCF and BAM files and generates a pileup for all the variants in the vcf file. See https://demichelislab.unitn.it/lib/exe/fetch.php?media=manual.pdf for ASEQ manual.
#' Two modes:
#'    -PILEUP: Default use. Generates a pileup from the BAM file
#'    -GENOTYPE: Generates genotype at VCF marked positions
#' The output is a ASEQ pileup format that can be used for downstream analysis. In addition, genotype mode generates a VCF file with all heterozygous calls.
#'
#'
#' @param bin_path [REQUIRED] Path to ASEQ binary. Default tools/ASEQ/binaries/linux64/ASEQ
#' @param vcf [REQUIRED] Path to vcf file.
#' @param bam [REQUIRED] Path to bam file.
#' @param mqr [OPTIONAL] Filter by read mapping quality >=. Default value 1. See ASEQ doc.
#' @param mbq [OPTIONAL] Filter by base quality >=. Default value 1. See ASEQ doc.
#' @param mdc [OPTIONAL] Filter by coverage per base >=. Default value 1. See ASEQ doc.
#' @param htperc [OPTIONAL] Heterozigosity test based on percentage. Default value 0.2. Only in GENOTYPE mode. See ASEQ doc.
#' @param pht [OPTIONAL] P-value for heterozigosity test. Default value 0.2. Only in GENOTYPE mode. See ASEQ doc.
#' @param mode [OPTIONAL] Usage mode. Default value PILEUP. See ASEQ doc.7
#' @param threads [OPTIONAL] Number of threads to use. Default value 3.
#' @param output_dir [OPTIONAL]  Path to the output directory.
#' @param verbose [OPTIONAL]  Enables progress messages. Default False.
#' @export



call_ASEQ=function(vcf="",bin_path="tools/ASEQ/binaries/linux64/ASEQ",bam="",mqr="",mbq="",mdc="",htperc="",pht="",mode="",output_dir="",threads=3,verbose=FALSE){

  sample_name=ULPwgs::get_sample_name(bam)

  sep="/"
  if(output_dir==""){
    sep=""
  }
  if (grepl(".gz",vcf)){
    system(paste("gunzip -d",vcf))
    vcf=sub(".gz","",vcf)
  }

  if (mode==""|mode=="PILEUP"){
    out_file_dir=paste0(output_dir,sep,sample_name,"_ASEQ_PILEUP")}
  else{
    out_file_dir=paste0(output_dir,sep,sample_name,"_ASEQ_GENOTYPING")}

  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  if (mode!=""){
    mode=paste0(" mode=",mode)
    if(mode=="GENOTYPE"){
        if (htperc!=""){
          htperc=paste0(" htperc=",htperc)
    }
        if (pht!=""){
          pht=paste0(" pht=",pht)
      }
    }
  }
  if (mqr!=""){
    mqr=paste0(" mqr=",mqr)
  }
  if (mbq!=""){
    mbq=paste0(" mbq=",mbq)
  }
  if (mdc!=""){
    mdc=paste0(" mdc=",mdc)
  }

  if(verbose){
    print(paste0(bin_path," vcf=",vcf, " bam=",bam,mode,mqr,mbq,mdc,htperc, " threads=",threads, " out=",out_file_dir))

  }
  system(paste0(bin_path," vcf=",vcf, " bam=",bam,mode,mqr,mbq,mdc,htperc, " threads=",threads, " out=",out_file_dir))

}

#' Format vcf for downstream point mutatation analysis.
#'
#'
#' This function takes a VCF file and outputs a tab delimited file with data for downstream analysis in CLONETv2 point mutations analysis.
#'
#' @param bin_path [REQUIRED] Path to ASEQ binary. Default tools/bcftools/bcftools
#' @param vcf_dir [REQUIRED] Path to vcf file directory.
#' @param pattern [OPTIONAL] Pattern to filter vcfs by.
#' @param vcf [OPTIONAL] Path to vcf file. Only if vcf directory not provided.
#' @param output_dir [OPTIONAL]  Path to the output directory.
#' @param output_name [OPTIONAL]  Name of output file.
#' @param verbose [OPTIONAL]  Enables progress messages. Default False.
#' @export


format_PM_data=function(bin_path="tools/bcftools/bcftools",vcf_dir="",pattern="",vcf="",output_dir="",output_name="",verbose=FALSE){




  if(vcf!="" & vcf_dir!=""){

    stop("Only vcf or vcf_dir arguments can be provided, not both.")
  }

  if (vcf_dir!=""){
    files=list.files(path=vcf_dir,pattern=pattern,full.names=TRUE)
  }else {
    files=vcf
  }

  if(output_name==""){
    output_name="PMreadsCount"
  }

  sep="/"
  if(output_dir==""){
    sep=""
  }

  out_file_dir=paste0(output_dir,sep,"FORMATTED_PM_ANALYSIS")
  out_file=paste0(out_file_dir,"/",output_name,".tsv")

  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  if(verbose){
    print(paste("printf","\'Gene.id\tchr\tstart\tend\tsample\tref.count\talt.count\n\'",">",out_file))
    lapply(files,FUN=function(x){print(paste0(bin_path," +", system(paste0("echo ",sub("bcftools$", "plugins\\1",bin_path),"/split-vep.so"),intern=TRUE)," -f \'%SYMBOL\\t%CHROM\\t%POS\\t[\\t%AD]\\n\' -i \'TYPE=\"snp\" & AD[0:1-]>1\' -s worst ",x," | awk -F \'[\t,]\' \'{ print  $1,$2,$3,$3, \"",ULPwgs::get_sample_name(x),"\", $5, $6 }\' OFS=\'\\t\' >> ",out_file))})
  }
    system(paste("printf","\'Gene.id\tchr\tstart\tend\tsample\tref.count\talt.count\n\'",">",out_file))
    lapply(files,FUN=function(x){system(paste0(bin_path," +",system(paste0("echo ",sub("bcftools$", "plugins\\1",bin_path),"/split-vep.so"),intern=TRUE), " -f \'%SYMBOL\\t%CHROM\\t%POS\\t[\\t%AD]\\n\' -i \'TYPE=\"snp\" & AD[0:1-]>1\' -s worst ",x," | awk -F \'[\t,]\' \'{ print  $1,$2,$3,$3, \"",ULPwgs::get_sample_name(x),"\", $5, $6 }\' OFS=\'\\t\' >> ",out_file))})
  }




#' Variant calling using Platypus
#'
#' This function calls somatic variants in a pair of tumor-normal matched samples, or
#' just in a tumor sample if no matched sample is not available.
#'
#' @param tumor_bam Path to tumor bam file.
#' @param normal_bam Path to germline bam file.
#' @param bin_path Path to fastQC executable. Default path tools/platypus/Platypus.py.
#' @param bin_path2 Path to bgzip binary. Default tools/htslib/bgzip.
#' @param bin_path3 Path to tabix binary. Default tools/htslib/tabix.
#' @param ref_genome Path to reference genome fasta file.
#' @param vcf_overlay Path to vcf overlay to use as source.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads to use. Default 3.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @export


call_platypus=function(bin_path="tools/platypus/Platypus.py",bin_path2="tools/htslib/bgzip",bin_path3="tools/htslib/tabix",tumor_bam="",normal_bam="",ref_genome="",vcf_overlay="",output_dir="",verbose=FALSE,threads=3,output_name=""){

  sep="/"

  if(output_dir==""){
    sep=""
  }

  if (output_name==""){
    sample_name=ULPwgs::get_sample_name(tumor_bam[1])
  }else{
    sample_name=output_name
  }


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
  out_file=paste0(out_file,"/",sample_name,".PLATYPUS.vcf")

# TO DO FIX THIS MESS

  if(verbose){
      print(paste0(bin_path," callVariants --refFile=",ref_genome, paste0(" --bamFiles=",tumor,",",norm)," --source=",vcf_overlay," --output=",out_file," --filterReadPairsWithSmallInserts=0 --minPosterior=0 --getVariantsFromBAMs=1 --logFileName=",paste0(out_file,".log")," --nCPU=",threads))

  }
  system(paste0(bin_path," callVariants --refFile=",ref_genome, paste0(" --bamFiles=",tumor,",",norm), " --source=",vcf_overlay," --output=",out_file," --filterReadPairsWithSmallInserts=0 --minPosterior=0 --getVariantsFromBAMs=1 --logFileName=",paste0(out_file,".log")," --nCPU=",threads))
  system(paste("cp", out_file, paste0(out_file,".tmp")))
  bgzip(bin_path=bin_path2,file=out_file)
  tab_indx(bin_path=bin_path3,file=paste0(out_file,".gz"))
  system(paste("cp", paste0(out_file,".tmp"), out_file))
  system(paste("rm -rf", paste0(out_file,".tmp")))
}
