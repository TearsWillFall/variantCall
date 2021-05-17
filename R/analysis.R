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


call_mutect2=function(region="",bin_path="tools/gatk/gatk",tumor_bam="",normal_bam="",ref_genome="",germ_resource="",pon="",output_dir="",output_name="",verbose=FALSE,orientation=TRUE){

  sep="/"

  if(output_dir==""){
    sep=""
  }

  if (output_name==""){
    sample_name=ULPwgs::get_sample_name(tumor_bam[1])
  }else{
    sample_name=output_name
  }

  out_file_dir=paste0(output_dir,sep,sample_name,"_MUTECT2_VARIANTS_VCF")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  reg=""
  if (region==""){
      out_file=paste0(out_file_dir,"/",sample_name,".UNFILTERED_MUTECT2.vcf")
  }else{
      reg=paste0(" -L ",region)
      out_file=paste0(out_file_dir,"/",sample_name,".",region,".UNFILTERED_MUTECT2.vcf")
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

  if (pon!=""){
    pon=paste0(" --panel-of-normals ",pon)
  }

  f1r2=""
  if (orientation){
      f1r2=paste0(" --f1r2-tar-gz ",out_file_dir,"/",sample_name,".",region,".f1r2.tar.gz")
  }

  if(verbose){
      print(paste0(bin_path," Mutect2 -R ",ref_genome,tumor, norm," --germline-resource ",germ_resource, pon," -O ",out_file, reg,f1r2))

  }
  system(paste0(bin_path," Mutect2 -R ",ref_genome,tumor, norm, " --germline-resource ",germ_resource, pon, " -O ",out_file, reg,f1r2))
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



#' Predict variant effect using VEP
#'
#' This function predicts the effect of variant found in a VCF file
#'
#' @param bin_path Path to bcftools binary. Default tools/ensembl-vep/vep.
#' @param bin_path2 Path to bgzip binary. Default tools/htslib/bgzip.
#' @param bin_path3 Path to bgzip binary. Default tools/htslib/tabix.
#' @param vcf Path to vcf file.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads to use. Default 3.
#' @export



call_vep=function(bin_path="tools/ensembl-vep/vep",bin_path2="tools/htslib/bgzip",bin_path3="tools/htslib/tabix",vcf="",verbose=FALSE,output_dir="",threads=3){
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

  out_file=paste0(out_file_dir,"/",sample_name,".VEP.vcf")

  if(verbose){
    print(paste(bin_path,"-i",vcf,"-o",out_file,"--cache --port 3337 --everything --force_overwrite --vcf --fork",threads))
  }
  system(paste(bin_path,"-i",vcf,"-o",out_file,"--cache --port 3337 --everything --force_overwrite --vcf --fork",threads))
  system(paste("cp", out_file, paste0(out_file,".tmp")))
  bgzip(bin_path=bin_path2,file=out_file)
  tab_indx(bin_path=bin_path3,file=paste0(out_file,".gz"))
  system(paste("cp", paste0(out_file,".tmp"), out_file))
  system(paste("rm -rf", paste0(out_file,".tmp")))
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
    print(paste(bin_path,"-i",vcf,"-o",out_file,"--cache --port 3337 --everything --force_overwrite --vcf --fork",threads))
  }
  system(paste(bin_path,"-i",vcf,"-o",out_file,"--cache --port 3337 --everything --force_overwrite --vcf --fork",threads))
}


#' Variant calling using MuTECT2 on parallel per genomic region
#'
#' This function calls somatic variants in a pair of tumor-normal matched samples, or
#' just in a tumor sample if no matched sample is not available.
#'
#' @param tumor_bam [OPTIONAL] Path to tumor bam file/s.
#' @param normal_bam [OPTIONAL] Path to germline bam file/s.
#' @param bam_dir [OPTIONAL] Path to bam_dir. Only if tumor_bam and/or normal_bam not given.
#' @param bin_path [REQUIRED] Path to gatk binary. Default path tools/gatk/gatk.
#' @param bin_path2 [REQUIRED] Path to bcftools binary. Default path tools/bcftools/bcftools.
#' @param bin_path3 [REQUIRED] Path to tabix binary. Default path tools/htslib/tabix.
#' @param bin_path4 [REQUIRED] Path to bgzip binary. Default path tools/htslib/bgzip.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param germ_resource [REQUIRED] Path to germline resources vcf file.
#' @param region_bed [REQUIRED] Path to bed file with regions to analyze.
#' @param patient_id [OPTIONAL] Patients ID to identify project with. Otherwise name of tumor sample will be used.
#' @param germ_pattern [OPTIONAL] Identifier for normal samples. Default GL. Only used when bam_dir is given.
#' @param pon [OPTIONAL] Path to panel of normal.
#' @param threads [OPTIONAL] Number of threads. Default 3
#' @param db [OPTIONAL] Path to vcf with known variants for contamination estimation.
#' @param interval [OPTIONAL] Path to vcf with intervals to analyze contamination.
#' @param chr_filter [OPTIONAL] Chromosomes to analyze. canonical/autosomal/all or a list of chromosomes
#' @param orientation [OPTIONAL] Generate a read orientation model to filter variants. Default False
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @export
#' @import pbapply


call_mutect2_parallel=function(bin_path="tools/gatk/gatk",bin_path2="tools/bcftools/bcftools",bin_path3="tools/htslib/bgzip",bin_path4="tools/htslib/tabix",tumor_bam="",normal_bam="",bam_dir="",patient_id="",germ_pattern="GL",ref_genome="",germ_resource="",pon="",output_dir="",region_bed="",db="",interval="",threads=3,verbose=FALSE,chr_filter="canonical",orientation=TRUE){

  if (patient_id==""){
    sample_name=ULPwgs::get_sample_name(tumor_bam[1])
  }else{
    sample_name=patient_id
  }

  sep="/"

  if(output_dir==""){
    sep=""
  }

  out_file_dir=paste0(output_dir,sep,sample_name,"_MUTECT2_VARIANTS_VCF")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  dat=read.table(region_bed)


  ## If a path to bam_dir is supplied and no individual bams are given identify bams in dir based on patient_id and germ_pattern

  if (tumor_bam=="" & normal_bam=="" & bam_dir!="" & patient_id!=""){
    files=list.files(bam_dir,recursive=TRUE,full.names=TRUE,pattern=patient_id)
    files=files[grepl("bam$",files)]
    tumor_bam=files[!grepl(germ_pattern,files)]
    normal_bam=files[grepl(germ_pattern,files)]
  }

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
  pbapply(X=dat[,c("Region"),drop=FALSE],1,FUN=call_mutect2,bin_path=bin_path,tumor_bam=tumor_bam,normal_bam=normal_bam,ref_genome=ref_genome,germ_resource=germ_resource,pon=pon,output_dir=output_dir,output_name=patient_id,verbose=verbose,orientation=orientation,cl=cl)
  on.exit(parallel::stopCluster(cl))

  vcf_concatenate(bin_path=bin_path2,vcf_dir=out_file_dir,output_dir=out_file_dir,verbose=verbose)
  vcf_sort(bin_path=bin_path2,vcf=paste0(out_file_dir,"/",sample_name,"_CONCATENATED","/",sample_name,".CONCATENATED.vcf"),output_dir=out_file_dir,verbose=verbose)
  vcf_stats_merge(bin_path=bin_path,vcf_stats_dir=out_file_dir,output_dir=out_file_dir,verbose=verbose)

  contamination_tables=""
  contamination_segments=""

  ## If common variant db and interval is given, then contamination will be estimated and included in the vcf filtering.

  if (db!="" & interval!=""){
    estimate_contamination_parallel(bin_path=bin_path,bam_dir=bam_dir,germ_pattern=germ_pattern,patient_id=patient_id,db=db,interval=interval,verbose=verbose,output_dir=out_file_dir,threads=threads)
    contamination_files=list.files(paste0(out_file_dir,"/",patient_id,"_CONTAMINATION"),full.names=TRUE)
    contamination_tables=contamination_files[grepl("contamination.table$",contamination_files)]
    contamination_segments=contamination_files[grepl("segments.table$",contamination_files)]
  }

  ## If orientation variable is established orientation_model will generated and used in vcf filtering

  if (orientation){
    learn_orientation(bin_path=bin_path,f1r2_dir=out_file_dir,output_name=sample_name,output_dir=out_file_dir,verbose=verbose)
  }

  vcf_filtering(bin_path=bin_path,bin_path2=bin_path3,bin_path3=bin_path4,ref_genome=ref_genome,unfil_vcf=paste0(out_file_dir,"/",patient_id,"_SORTED/",patient_id,".SORTED.CONCATENATED.vcf"),
  unfil_vcf_stats=paste0(out_file_dir,"/",patient_id,"_MERGED_VCF_STATS/",sample_name,".MERGED.vcf.stats"),output_dir=out_file_dir,verbose=verbose,contamination=contamination_tables,segmentation=contamination_segments,
  orientation_model=paste0(out_file_dir,"/",patient_id,"_ORIENTATION_MODEL/",sample_name,".read-orientation-model.tar.gz"))

  ## Remove intermediary files
  system(paste0("rm ",out_file_dir,"/*"))
  system(paste0("rm -rf ",out_file_dir,"/",sample_name,"_CONCATENATED"))
  system(paste0("rm -rf ",out_file_dir,"/",sample_name,"_MERGED_VCF_STATS"))
  system(paste0("rm -rf ",out_file_dir,"/",sample_name,"_SORTED"))
}

#' Wrapper for GATK HaplotypeCaller for Germline Variant calling
#'
#' This function takes an annotated VCF with CNN_D1 or CNN_D2 scores and filter the
#' variants based on set threshold.Cuttently only single sample implementation
#'
#' @param bin_path [REQUIRED] Path to GATK binary. Default tools/gatk/gatk
#' @param normal_bam [REQUIRED] Path to BAM file
#' @param ref_genome [REQUIRED] Path to reference genome
#' @param resources [OPTIONAL] Path to resources for variant filtering
#' @param region [OPTIONAL] Genomic region/regions to analyze. In samtools format.
#' @param output_name [OPTIONAL] Name of the sample to output
#' @param info_key [OPTIONAL] Annotation column to select. Default CNN_D1
#' @param snp_tranche [OPTIONAL] SNP tranche filter value. Default 99.95
#' @param indel_tranche [OPTIONAL] Indel tranche filter value. Default 99.4
#' @param keep_previous_filters [OPTIONAL] Keep previous filters in VCF. Default False
#' @param output_dir [OPTIONAL] Path to output dir
#' @param patient_id [OPTIONAL] Patient ID after which to name the output. If not given bam file name will be used.
#' @param threads [OPTIONAL] Number of threads per job. Default 3
#' @param verbose [Optional] Enables progress messages. Default False
#' @export

call_HaplotypeCaller=function(bin_path="tools/gatk/gatk",normal_bam="",ref_genome="",region="",output_dir="",resources="",patient_id="",info_key="CNN_1D",snp_tranche=99.95,indel_tranche=99.4,keep_previous_filters=FALSE,verbose=FALSE,threads=3){

  sep="/"
  if(output_dir==""){
    sep=""
  }

  if (patient_id==""){
    sample_name=ULPwgs::get_sample_name(normal_bam[1])
  }else{
    sample_name=patient_id
  }

  out_file_dir=paste0(output_dir,sep,sample_name,"_HAPLOTYPECALLER_VARIANTS_VCF")

  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  if(region!=""){
    region=paste0(" -L ",region,collapse=" -L ")
  }

  if(verbose){
    print(paste(bin_path," HaplotypeCaller -R ",ref_genome," -I ",normal_bam," -O ",paste0(out_file_dir,"/",sample_name,".GL.vcf.gz ")," --native-pair-hmm-threads ",threads,region))
  }
  system(paste(bin_path," HaplotypeCaller -R ",ref_genome," -I ",normal_bam," -O ",paste0(out_file_dir,"/",sample_name,".GL.vcf.gz ")," --native-pair-hmm-threads ",threads,region))


  if (info_key=="CNN_1D"){
      CNNScoreVariants(bin_path=bin_path,vcf=paste0(out_file_dir,"/",sample_name,".GL.vcf.gz "),ref_genome=ref_genome,output_dir=out_file_dir,output_name=sample_name,verbose=verbose)
      scored_vcf=paste0(out_file_dir,"/",sample_name,"_CNNscored","/",sample_name,".CNNscored.1D.vcf")
  }else if(info_key=="CNN_2D"){
      CNNScoreVariants(bin_path=bin_path,vcf=paste0(out_file_dir,"/",sample_name,".GL.vcf.gz "),ref_genome=ref_genome,bam=normal_bam,output_dir=out_file_dir,output_name=sample_name,verbose=verbose)
      scored_vcf=paste0(out_file_dir,"/",sample_name,"_CNNscored","/",sample_name,".CNNscored.2D.vcf")
  }

  FilterVariantTranches(bin_path=bin_path,vcf=scored_vcf,resources=resources,output_name=sample_name,info_key=info_key,
  snp_tranche=snp_tranche,indel_tranche=indel_tranche,output_dir=out_file_dir,keep_previous_filters=keep_previous_filters,verbose=verbose)
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






#' Variant calling of Somatic and Germline (MuTECT2 and Platypus) mutations and produces their subsequent annotation (VEP)
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
#' @param db [OPTIONAL] Path to vcf with common variants. Used for contamination estimation.
#' @param interval [OPTIONAL] Path to interval for common variants to analyze. Used for contamination estimation.
#' @param orientation [OPTIONAL] Generate a read orientation model to filter variants. Default TRUE
#' @param info_key [OPTIONAL] Annotation column to select. Default CNN_D2
#' @param snp_tranche [OPTIONAL] SNP tranche filter value. Default 99.95
#' @param indel_tranche [OPTIONAL] Indel tranche filter value. Default 99.44
#' @param targeted [OPTIONAL] Sequencing method Exome/Targeted. Default TRUE
#' @param chr_filter [OPTIONAL] Chromosomes to analyze. canonical/autosomal/all or a list of chromosomes
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @export

call_variants=function(bin_path="tools/gatk/gatk",bin_path2="tools/bcftools/bcftools",bin_path3="tools/htslib/bgzip",bin_path4="tools/htslib/tabix",bin_path5="tools/platypus/Platypus.py",bin_path6="tools/ensembl-vep/vep",bam_dir="",patient_id="",germ_pattern="GL",ref_genome="",
germ_resource="",pon="",output_dir="",region_bed="",chr_filter="canonical",db="",interval="",orientation=TRUE,resources="",info_key="CNN_2D",snp_tranche=99.95,indel_tranche=99.4,threads=3,verbose=FALSE,targeted=TRUE){

    sep="/"
    if(output_dir==""){
      sep=""
    }
    out_file_dir=paste0(output_dir,sep,patient_id,"_VARIANTS")
    if (!dir.exists(out_file_dir)){
        dir.create(out_file_dir)
    }

    chr_pass=""
    if(chr_filter=="canonical"){
      chr_pass=c(1:22,"X","Y","M","MT")
    ## Include autosomal chromosomal and nothing else
    }else if(chr_filter=="autosomal"){
      chr_pass=c(1:22)
    ## Include only specific chromosomes
    }else if(chr_filter!="all"){
      chr_pass=chr_filter
    }

    files=list.files(bam_dir,recursive=TRUE,full.names=TRUE,pattern=patient_id)
    files=files[grepl("bam$",files)]
    tumor_bam=files[!grepl(germ_pattern,files)]
    normal_bam=files[grepl(germ_pattern,files)]
    ## Call Somatic SNVs+INDELs Using Mutect2
    call_mutect2_parallel(bin_path=bin_path,bin_path2=bin_path2,bin_path3=bin_path3,bin_path4=bin_path4,tumor_bam=tumor_bam,normal_bam=normal_bam,bam_dir=bam_dir,germ_pattern=germ_pattern,ref_genome=ref_genome,germ_resource=germ_resource,pon=pon,output_dir=out_file_dir,region_bed=region_bed,threads=threads,verbose=verbose,patient_id=patient_id,chr_filter=chr_filter,orientation=orientation,interval=interval,db=db)
    ## Call Germline SNVs+INDELs Using HaplotypeCaller
    call_HaplotypeCaller(bin_path=bin_path,normal_bam=normal_bam,ref_genome=ref_genome,output_dir=out_file_dir,resources=resources,info_key=info_key,snp_tranche=snp_tranche,indel_tranche=indel_tranche,patient_id=patient_id,verbose=verbose,threads=threads,region=chr_pass)
    ## Call Germline + Somatic SNVs+INDELs Using Platypus
    call_platypus(bin_path=bin_path5,bin_path2=bin_path3,bin_path3=bin_path4,tumor_bam=tumor_bam,normal_bam=normal_bam,ref_genome=ref_genome,vcf_overlay=paste0(out_file_dir,"/",patient_id,"_MUTECT2_VARIANTS_VCF/",patient_id,"_FILTERED/",patient_id,".FILTERED.vcf.gz"),output_dir=out_file_dir,verbose=verbose,threads=threads,output_name=patient_id,targeted=targeted)
    ## Anotate Variants using VEP
    call_vep(bin_path=bin_path6,bin_path2=bin_path3,bin_path3=bin_path4,vcf=paste0(out_file_dir,"/",patient_id,"_HAPLOTYPECALLER_VARIANTS_VCF/",patient_id,"_FILTERED_TRENCHES/",patient_id,".FILTERED.vcf.gz"),verbose=verbose,output_dir=paste0(out_file_dir,"/",patient_id,"_HAPLOTYPECALLER_VARIANTS_VCF"),threads=threads)
    call_vep(bin_path=bin_path6,bin_path2=bin_path3,bin_path3=bin_path4,vcf=paste0(out_file_dir,"/",patient_id,"_MUTECT2_VARIANTS_VCF/",patient_id,"_FILTERED/",patient_id,".FILTERED.vcf.gz"),verbose=verbose,output_dir=paste0(out_file_dir,"/",patient_id,"_MUTECT2_VARIANTS_VCF"),threads=threads)
    call_vep(bin_path=bin_path6,bin_path2=bin_path3,bin_path3=bin_path4,vcf=paste0(out_file_dir,"/",patient_id,"_PLATYPUS_VARIANTS_VCF/",patient_id,".PLATYPUS.vcf.gz"),verbose=verbose,output_dir=paste0(out_file_dir,"/",patient_id,"_PLATYPUS_VARIANTS_VCF"),threads=threads)
}


#' Filters SNVs from any variant caller to remove false positives (FP) using FiNGS on parallel
#'
#' This function takes a path to a directory with BAM files and a patient ID. Then it subsets all BAM
#' files specific to the patient, identifying between cancer and normal samples using the germline identifier.
#' It also takes a directory with vcf files and subsets them based on the sample ID.
#' Afterwards, it calls FiNGS to calculates metrics based on BAM files, providing filtering of FP.
#'
#'
#' @param bin_path [REQUIRED] Path to FiNGS binary. Default tools/FiNGS/fings/FiNGS.py.
#' @param bin_path2 [REQUIRED] Path to bcftools binary. Default tools/bcftools/bcftools.
#' @param bin_path3 [REQUIRED] Path to bgzip binary. Default tools/htslib/bgzip.
#' @param bin_path4 [REQUIRED] Path to tabix binary. Default tools/htslib/tabix.
#' @param bam_dir [REQUIRED] Path to directory with BAM files.
#' @param vcf_dir [REQUIRED] Path to directory with VCF files.
#' @param patient_id [REQUIRED] Patient ID to analyze. Has to be in file names to subselect samples.
#' @param germ_pattern [REQUIRED] Pattern used to identify germline samples. Ex GL
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param max_depth [OPTIONAL] Maximum number of reads.Reads beyond this depth will be ignored
#' @param pass_in [OPTIONAL] Only input variants with a PASS. Default TRUE.
#' @param pass_out [OPTIONAL] Only output variants with a PASS. Default FALSE.
#' @param param [REQUIRED] Path to file with filter params.
#' @param threads [OPTIONAL] Number of threads. Default 3
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @export

call_fings_parallel=function(bin_path="tools/FiNGS/fings/FiNGS.py",bin_path2="tools/bcftools/bcftools",bin_path3="tools/htslib/bgzip",bin_path4="tools/htslib/tabix",bam_dir="",vcf_dir="",patient_id="",germ_pattern="GL",ref_genome="",max_depth=1000,pass_in=TRUE,pass_out=FALSE,param="tools/fings/icgc_filter_parameters.txt",threads=3,output_dir="",verbose=FALSE){
    sep="/"
    if(output_dir==""){
      sep=""
    }
    out_file_dir=paste0(output_dir,sep,patient_id,"_FiNGS")
    if (!dir.exists(out_file_dir)){
        dir.create(out_file_dir)
    }
    files=list.files(bam_dir,recursive=TRUE,full.names=TRUE,pattern=patient_id)
    files=files[grepl("bam$",files)]
    tumor_bams=files[!grepl(germ_pattern,files)]
    normal_bam=files[grepl(germ_pattern,files)]

    files=list.files(vcf_dir,recursive=TRUE,full.names=TRUE,pattern=patient_id)
    files=files[grepl("vcf$",files)]
    tumor_vcfs=files[!grepl(germ_pattern,files)]
    tumor_bams=data.frame(bam=tumor_bams,sample_id=sapply(tumor_bams,FUN=ULPwgs::get_sample_name),stringsAsFactors =FALSE)
    tumor_vcfs=data.frame(vcf=tumor_vcfs,sample_id=sapply(tumor_vcfs,FUN=ULPwgs::get_sample_name),stringsAsFactors =FALSE)
    files=dplyr::left_join(tumor_bams,tumor_vcfs,by="sample_id")
    if (any(is.na(files))){
      stop("Could not match BAM and VCF file IDs.")
    }
    cl=parallel::makeCluster(threads)
    pbapply::pblapply(X=1:nrow(files),FUN=function(x){call_fings(bin_path=bin_path,bin_path2=bin_path2,bin_path3=bin_path3,bin_path4=bin_path4,tumor_bam=files[x,]$bam,normal_bam=normal_bam,tumor_vcf=files[x,]$vcf,ref_genome=ref_genome,output_dir=out_file_dir,max_depth=max_depth,verbose=verbose,param=param)},cl=cl)
    on.exit(parallel::stopCluster(cl))

    files=list.files(out_file_dir,recursive=TRUE,full.names=TRUE,pattern="tumor.combined.txt")
    cl=parallel::makeCluster(threads)
    dat=pbapply::pblapply(X=files,FUN=function(x){dat=read.table(x,header=TRUE);dat$sample=ULPwgs::get_sample_name(sub("_FiNGS_FILTERED/tumor.combined.txt","",x));return(dat)},cl=cl)
    on.exit(parallel::stopCluster(cl))
    dat=dat %>% dplyr::bind_rows()
    write.table(dat,paste0(out_file_dir,"/",patient_id,".tumor.stats.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE)
    files=list.files(out_file_dir,recursive=TRUE,full.names=TRUE,pattern="tumor.combined.txt")

    files=list.files(out_file_dir,recursive=TRUE,full.names=TRUE,pattern="normal.combined.txt")
    cl=parallel::makeCluster(threads)
    dat=pbapply::pblapply(X=files,FUN=function(x){dat=read.table(x,header=TRUE);dat$sample=ULPwgs::get_sample_name(sub("_FiNGS_FILTERED/normal.combined.txt","",x));return(dat)},cl=cl)
    on.exit(parallel::stopCluster(cl))
    dat=dat %>% dplyr::bind_rows()
    write.table(dat,paste0(out_file_dir,"/",patient_id,".normal.stats.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE)
  }

#' Filters SNVs from any variant caller to remove false positives (FP) using FiNGS
#'
#' This function takes a paths of a pair of matched tumor/normal BAMs and the tumors' VCF.
#' Then it calls FiNGS to calculate metrics based on the BAM files, as well as providing filtering
#' for FP variants called by the variant caller.
#'
#'
#' @param bin_path [REQUIRED] Path to FiNGS binary. Default tools/FiNGS/fings/FiNGS.py.
#' @param bin_path2 [REQUIRED] Path to bcftools binary. Default tools/bcftools/bcftools.
#' @param bin_path3 [REQUIRED] Path to bgzip binary. Default tools/htslib/bgzip.
#' @param bin_path4 [REQUIRED] Path to tabix binary. Default tools/htslib/tabix.
#' @param bam_dir [REQUIRED] Path to directory with BAM files.
#' @param vcf_dir [REQUIRED] Path to directory with VCF files.
#' @param patient_id [REQUIRED] Patient ID to analyze. Has to be in file names to subselect samples.
#' @param germ_pattern [REQUIRED] Pattern used to identify germline samples. Ex GL
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param max_depth [OPTIONAL] Maximum number of reads.Reads beyond this depth will be ignored
#' @param pass_in [OPTIONAL] Only input variants with a PASS. Default TRUE.
#' @param pass_out [OPTIONAL] Only output variants with a PASS. Default FALSE.
#' @param param [REQUIRED] Path to file with filter params.
#' @param threads [OPTIONAL] Number of threads. Default 3
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @export


call_fings=function(bin_path="tools/fings/FiNGS.py",bin_path2="tools/bcftools/bcftools",bin_path3="tools/htslib/bgzip",bin_path4="tools/htslib/tabix",tumor_bam="",normal_bam="",tumor_vcf="",ref_genome="",output_dir="",max_depth=1000,pass_in=TRUE,pass_out=FALSE,verbose=FALSE,param="tools/fings/icgc_filter_parameters.txt"){

  sep="/"
  if(output_dir==""){
    sep=""
  }

  sample_name=ULPwgs::get_sample_name(tumor_bam)
  out_file_dir=paste0(output_dir,sep,sample_name,"_FiNGS_FILTERED")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  pass_i=""
  if (pass_in){
    pass_i=" --PASSonlyin "
  }
  pass_o=""
  if (pass_out){
    pass_o=" --PASSonlyout "
  }

  if(verbose){
    print(paste(bin_path," -n ",normal_bam," -t ", tumor_bam," -v ", tumor_vcf, " -m ",max_depth, " -r ",ref_genome, " -p ",param,pass_i,pass_o, " -d ", out_file_dir))
  }
  system(paste(bin_path," -n ",normal_bam," -t ", tumor_bam," -v ", tumor_vcf, " -m ",max_depth, " -r ",ref_genome, " -p ",param,pass_i,pass_o," -d ",  out_file_dir))
  vcf_filter_variants(unfil_vcf=paste0(out_file_dir,"/",ULPwgs::get_sample_name(tumor_bam),".filtered.vcf"),bin_path=bin_path2,bin_path2=bin_path3,bin_path3=bin_path4,qual="",mq="",state="",ref="",type="",filter="PASS",verbose=verbose,output_dir=out_file_dir)
  system(paste0("gunzip ",out_file_dir,"/*.gz"))
}







#' Call Somatic/Germline single nucleotide variants using Strelka

#' This function takes a pair of matched samples or single | multiple germline variants
#' and calls variants on them. To specify the variant calling mode supply the right strelka binary.
#' It is recommended to supply a vcf with indel candidates. This can be generated using MANTA workflow
#'
#' @param bin_path [REQUIRED] Path to strelka binary. Somatic or Germline.
#' @param tumor_bam [REQUIRED] Path to tumor  bam file.
#' @param normal_bam [OPTIONAL] Path to normal samples bam files.
#' @param ref_genome [REQUIRED] Path to reference genome.
#' @param indel_candidates [OPTIONAL] Path to indel candidates produced by MANTA.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param targeted [OPTIONAL] If exome/capture method. Default FALSE
#' @param threads [OPTIONAL] Number of threads per job. Default 3
#' @param exec_options [OPTIONAL] Type of execution. local (Single node) / sge (multiple nodes). Default local.
#' @param verbose [DEFAULT==FALSE] Enables progress messages.
#' @export

call_variants_strelka=function(bin_path="tools/strelka-2.9.10/build/bin/configureStrelkaSomaticWorkflow.py",tumor_bam="",normal_bam="",ref_genome="",indel_candidates="",output_dir="",verbose=FALSE,targeted=FALSE,threads=3,exec_options="local"){

  sep="/"
  if(output_dir==""){
    sep=""
  }

  if (tumor_bam!=""){
    sample_name=ULPwgs::get_sample_name(tumor_bam)
    out_file_dir=paste0(output_dir,sep,sample_name,"_STRELKA_SNV_SOMATIC")
    tumor_bam=paste(" --tumorBam ",tumor_bam)
    normal_bam=paste(" --normalBam ",normal_bam)
  }else{
    sample_name=ULPwgs::get_sample_name(normal_bam)
    out_file_dir=paste0(output_dir,sep,sample_name,"_STRELKA_SNV_GERMLINE")
    normal_bam=paste0(" --bam ", normal_bam,collapse=" --bam ")
  }

  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  if (indel_candidates!=""){
    indel_candidates=paste(" --indelCandidates ",indel_candidates)
  }

  exome=""
  if (targeted){
    exome=" --exome "
  }



  if(verbose){
    print(paste("python2.7 ",bin_path, tumor_bam, normal_bam, " --referenceFasta ", ref_genome,indel_candidates," --runDir ",out_file_dir, exome))
  }
  system(paste("python2.7 ",bin_path, tumor_bam, normal_bam, " --referenceFasta ", ref_genome,indel_candidates," --runDir ",out_file_dir, exome))


  if(verbose){
    print(paste0("python2.7 ",out_file_dir,"/runWorkflow.py -m ",exec_options," -j ",threads))
  }
  system(paste0("python2.7 ",out_file_dir,"/runWorkflow.py -m ",exec_options," -j ",threads))
}



#' Call Somatic/Germline structural variants using Manta

#' This function takes a pair of matched samples or single | multiple germline variants
#' and calls variants on them. Variant calling mode is established based on wether tumor_bam is provided or not.
#'
#' @param bin_path [REQUIRED] Path to strelka binary. Somatic or Germline.
#' @param tumor_bam [REQUIRED] Path to tumor  bam file.
#' @param normal_bam [OPTIONAL] Path to normal samples bam files.
#' @param ref_genome [REQUIRED] Path to reference genome.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param targeted [OPTIONAL] If exome/capture method. Default FALSE
#' @param threads [OPTIONAL] Number of threads per job. Default 3
#' @param verbose [DEFAULT==FALSE] Enables progress messages.
#' @export


call_sv_manta=function(bin_path="tools/manta-1.6.0/build/bin/configManta.py",tumor_bam="",normal_bam="",ref_genome="",output_dir="",verbose=FALSE,targeted=FALSE,threads=3){

  sep="/"
  if(output_dir==""){
    sep=""
  }

  if (tumor_bam!=""){
    sample_name=ULPwgs::get_sample_name(tumor_bam)
    tumor_bam=paste(" --tumorBam ",tumor_bam)
    normal_bam=paste(" --normalBam ",normal_bam)
    out_file_dir=paste0(output_dir,sep,sample_name,"_MANTA_SV_SOMATIC")
  }else{
    sample_name=ULPwgs::get_sample_name(normal_bam)
    normal_bam=paste0(" --bam ", normal_bam,collapse=" --bam ")
    out_file_dir=paste0(output_dir,sep,sample_name,"_MANTA_SV_GERMLINE")
  }

  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  exome=""
  if (targeted){
    exome=" --exome "
  }

  if(verbose){
    print(paste("python2.7 ",bin_path, tumor_bam, normal_bam, " --referenceFasta ", ref_genome," --runDir ",out_file_dir, exome))
  }
  system(paste("python2.7 ",bin_path, tumor_bam, normal_bam, " --referenceFasta ", ref_genome," --runDir ",out_file_dir, exome))


  if(verbose){
    print(paste0("python2.7 ",out_file_dir,"/runWorkflow.py -j ",threads))
  }
  system(paste0("python2.7 ",out_file_dir,"/runWorkflow.py  -j ",threads))
}



#' Call Somatic/Germline variants using Strelka in parallel

#' This function takes a pair of matched samples or single | multiple germline variants
#' and calls variants on them. To specify the variant calling mode supply the right strelka binary.
#' It is recommended to supply a vcf with indel candidates. This can be generated using MANTA workflow.
#' This function parallelizes the number of jobs that can be run.
#'
#' @param bin_path [REQUIRED] Path to strelka binary. Somatic or Germline.
#' @param bam_dir [REQUIRED] Path to directory with BAM files.
#' @param ref_genome [REQUIRED] Path to reference genome.
#' @param patient_id [REQUIRED] Patient ID to analyze. Has to be in file names to subselect samples.
#' @param germ_pattern [REQUIRED] Pattern used to identify germline samples. Ex GL
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param targeted [OPTIONAL] If exome/capture method. Default FALSE
#' @param threads [OPTIONAL] Number of threads per job. Default 3
#' @param jobs [OPTIONAL] Number of jobs. Default 1
#' @param exec_options [OPTIONAL] Type of execution. local (Single node) / sge (multiple nodes). Default local.
#' @param verbose [DEFAULT==FALSE] Enables progress messages.
#' @export


call_variants_strelka_parallel=function(bin_path="tools/strelka-2.9.10/build/bin/configureStrelkaSomaticWorkflow.py",bin_path2="~/tools/manta-1.6.0/build/bin/configManta.py",bam_dir="",ref_genome="",output_dir="",patient_id="",germ_pattern="GL",verbose=FALSE,targeted=FALSE,jobs=1,threads=3,exec_options="local"){

  sep="/"
  if(output_dir==""){
    sep=""
  }

  out_file_dir=paste0(output_dir,sep,patient_id,"_STRELKA_VARIANTS_VCF")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  files=list.files(bam_dir,recursive=TRUE,full.names=TRUE,pattern=patient_id)
  files=files[grepl("bam$",files)]
  tumor_bams=files[!grepl(germ_pattern,files)]
  normal_bam=files[grepl(germ_pattern,files)]

  cl=parallel::makeCluster(jobs)
  pbapply::pblapply(X=1:length(tumor_bams),FUN=function(x){
    call_sv_manta(bin_path=bin_path2,tumor_bam=tumor_bams[x],normal_bam=normal_bam,ref_genome=ref_genome,output_dir=out_file_dir,verbose=verbose,targeted=targeted,threads=threads);
    call_variants_strelka(bin_path=bin_path,tumor_bam=tumor_bams[x],normal_bam=normal_bam,ref_genome=ref_genome,output_dir=out_file_dir,verbose=verbose,indel_candidates=paste0(out_file_dir,"/",ULPwgs::get_sample_name(x),"_MANTA_SV_SOMATIC/results/variants/candidateSmallIndels.vcf.gz"),targeted=targeted,threads=threads,exec_options=exec_options)},cl=cl)
  on.exit(parallel::stopCluster(cl))
}


#' Call Somatic/Germline variants using MANTA in parallel

#' This function takes a pair of matched samples or single | multiple germline variants
#' and calls variants on them.
#' This function parallelizes the number of jobs that can be run at single time.
#'
#' @param bin_path [REQUIRED] Path to strelka binary. Somatic or Germline.
#' @param bam_dir [REQUIRED] Path to directory with BAM files.
#' @param ref_genome [REQUIRED] Path to reference genome.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param targeted [OPTIONAL] If exome/capture method. Default FALSE
#' @param threads [OPTIONAL] Number of threads per job. Default 3
#' @param jobs [OPTIONAL] Number of jobs. Default 1
#' @param verbose [DEFAULT==FALSE] Enables progress messages.
#' @export


call_sv_manta_parallel=function(bin_path="tools/manta-1.6.0/build/bin/configManta.py",bam_dir="",ref_genome="",output_dir="",verbose=FALSE,targeted=FALSE,jobs=1,threads=3){

  sep="/"
  if(output_dir==""){
    sep=""
  }

  out_file_dir=paste0(output_dir,sep,patient_id,"_MANTA_VARIANTS_VCF")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  files=list.files(bam_dir,recursive=TRUE,full.names=TRUE,pattern=patient_id)
  files=files[grepl("bam$",files)]
  tumor_bams=files[!grepl(germ_pattern,files)]
  normal_bam=files[grepl(germ_pattern,files)]

  cl=parallel::makeCluster(jobs)
  pbapply::pblapply(X=1:length(tumor_bams),FUN=function(x){call_sv_manta(bin_path=bin_path,tumor_bam=tumor_bams[x],normal_bam=normal_bam,ref_genome=ref_genome,output_dir=out_file_dir,verbose=verbose,targeted=targeted,threads=threads)},cl=cl)
  on.exit(parallel::stopCluster(cl))
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
#' @param mrq [OPTIONAL] Filter by read mapping quality >=. Default value 1. See ASEQ doc.
#' @param mbq [OPTIONAL] Filter by base quality >=. Default value 1. See ASEQ doc.
#' @param mdc [OPTIONAL] Filter by coverage per base >=. Default value 1. See ASEQ doc.
#' @param htperc [OPTIONAL] Heterozigosity test based on percentage. Default value 0.2. Only in GENOTYPE mode. See ASEQ doc.
#' @param pht [OPTIONAL] P-value for heterozigosity test. Default value 0.2. Only in GENOTYPE mode. See ASEQ doc.
#' @param mode [OPTIONAL] Usage mode. Default value PILEUP. See ASEQ doc.7
#' @param threads [OPTIONAL] Number of threads to use. Default value 3.
#' @param output_dir [OPTIONAL]  Path to the output directory.
#' @param verbose [OPTIONAL]  Enables progress messages. Default False.
#' @export



call_ASEQ=function(vcf="",bin_path="tools/ASEQ/binaries/linux64/ASEQ",bam="",mrq="",mbq="",mdc="",htperc="",pht="",mode="",output_dir="",threads=3,verbose=FALSE){

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
  if (mrq!=""){
    mrq=paste0(" mrq=",mrq)
  }
  if (mbq!=""){
    mbq=paste0(" mbq=",mbq)
  }
  if (mdc!=""){
    mdc=paste0(" mdc=",mdc)
  }

  if(verbose){
    print(paste0(bin_path," vcf=",vcf, " bam=",bam,mode,mrq,mbq,mdc,htperc, " threads=",threads, " out=",out_file_dir))

  }
  system(paste0(bin_path," vcf=",vcf, " bam=",bam,mode,mrq,mbq,mdc,htperc, " threads=",threads, " out=",out_file_dir))

}

#' Structural variant calling using svaba
#'
#' This function calls structural variants in a pair of tumor-normal matched samples
#' using svaba
#'
#'
#' @param tumor_bam [REQUIRED]  Path to tumor bam file.
#' @param normal_bam [OPTIONAL] Path to germline bam file.
#' @param bin_path [REQUIRED] Path to svaba binary executable. Default path tools/svaba/svaba.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param threads [OPTIONAL] Number of threads to use. Default 3.
#' @param targets [OPTIONAL] BED file with capture target regions.
#' @param dbsnp_indels [OPTIONAL] Database with indel annotations.
#' @param output_name [OPTIONAL] Name for the output.
#' @export

call_sv_svaba=function(tumor_bam="",bin_path="tools/svaba/bin/svaba",normal_bam="",ref_genome="",threads=3,output_name="",targets="",dbsnp_indels="",verbose=FALSE,output_dir=""){
  sep="/"

  if(output_dir==""){
    sep=""
  }

  if (output_name==""){
    sample_name=ULPwgs::get_sample_name(tumor_bam[1])
  }else{
    sample_name=output_name
  }

  out_file_dir=paste0(output_dir,sep,sample_name,"_SV_SVABA")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  out_file=paste0(out_file_dir,"/",sample_name)

  tgs=""
  if (!targets==""){
    tgs=paste0(" -k ",targets)
  }
  if(length(tumor_bam)>1){
    tumor_bam=paste0(tumor_bam,collapse=" -t ")
  }

  norm=""
  if (!normal_bam==""){
    if (length(normal_bam)>1){
        norm=paste0(" -n ",paste(normal_bam,collapse=" -n "))
    }else{
        norm=paste0(" -n ",normal_bam)
    }

  }

  dbsnp=""
  if (!dbsnp_indels==""){
    dbsnp=paste0(" -D ",dbsnp_indels)
  }

  if(verbose){
      print(paste0(bin_path," run -t ",tumor_bam,norm,tgs," -a ",out_file," -p ",threads," -G ",ref_genome,dbsnp))
  }
    system(paste0(bin_path," run -t ",tumor_bam,norm,tgs," -a ",out_file," -p ",threads," -G ",ref_genome,dbsnp))
}


#' Structural variant calling using svaba in parallel
#'
#' This function calls structural variants in a pair of tumor-normal matched samples pairs
#' using svaba in parallel
#'
#'
#' @param bam_dir [REQUIRED]  Path to bam files.
#' @param bin_path [REQUIRED] Path to svaba binary executable. Default path tools/svaba/svaba.
#' @param patient_id [REQUIRED] Patient ID patter to filter samples with.
#' @param germ_pattern [OPTIONAL] Patient ID pattern to filter samples with.
#' @param ref_genome [REQUIRED] Germline pattern to filter samples with
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param jobs [OPTIONAL] Number of paired jobs to run. Only when par_type is Paired. Default 3.
#' @param threads [OPTIONAL] Number of threads to use per job. Default 1.
#' @param targets [OPTIONAL] BED file with capture target regions.
#' @param par_type [OPTIONAL] Parallelization type. Joint: Runs a single job with multiple samples jointly | Paired: Runs a pair of tumor-normal through multiple jobs.
#' @export


call_sv_svaba_parallel=function(bin_path="tools/svaba/bin/svaba",targets="",bam_dir="",patient_id="",germ_pattern="GL",ref_genome="",jobs=3,threads=1,par_type="Joint",verbose=FALSE,output_dir=""){
  sep="/"
  if(output_dir==""){
    sep=""
  }
  out_file_dir=paste0(output_dir,sep,patient_id,"_SVABA_VARIANTS_VCF")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  files=list.files(bam_dir,recursive=TRUE,full.names=TRUE,pattern=patient_id)
  files=files[grepl("bam$",files)]
  tumor_bams=files[!grepl(germ_pattern,files)]
  normal_bam=files[grepl(germ_pattern,files)]

  if(par_type=="Joint"){
    call_sv(tumor_bam=tumor_bams,bin_path=bin_path,normal_bam=normal_bam,ref_genome=ref_genome,threads=threads,verbose=verbose,output_dir=out_file_dir,output_name=patient_id)
  }else if(par_type=="Paired"){
    cl=parallel::makeCluster(jobs)
    dat=pbapply::pblapply(X=tumor_bams,FUN=call_sv,bin_path=bin_path,normal_bam=normal_bam,ref_genome=ref_genome,threads=threads,cl=cl,verbose=verbose,output_dir=out_file_dir)
    on.exit(parallel::stopCluster(cl))
  }
}


#' Variant calling using Platypus
#'
#' This function calls somatic variants in a pair of tumor-normal matched samples, or
#' just in a tumor sample if no matched sample is available.
#'
#' @param tumor_bam [REQUIRED] Path to tumor bam file.
#' @param normal_bam [REQUIRED] Path to germline bam file.
#' @param bin_path [REQUIRED] Path to fastQC executable. Default path tools/platypus/Platypus.py
#' @param bin_path2 [REQUIRED] Path to bgzip binary. Default tools/htslib/bgzip.
#' @param bin_path3 [REQUIRED] Path to tabix binary. Default tools/htslib/tabix.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param vcf_overlay [REQUIRED] Path to vcf overlay to use as source.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL]  Enables progress messages. Default False.
#' @param threads [OPTIONAL]  Number of threads to use. Default 3.
#' @param targeted [OPTIONAL]  Sequencing data is exome/targeted. Default FALSE
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first sample in alphanumerical order will be used.
#' @export

call_platypus=function(bin_path="tools/platypus/Platypus.py",bin_path2="tools/htslib/bgzip",bin_path3="tools/htslib/tabix",tumor_bam="",normal_bam="",ref_genome="",vcf_overlay="",output_dir="",verbose=FALSE,threads=3,output_name="",targeted=FALSE){

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


  if (targeted){
    opt=paste0(" --filterDuplicates=0 ")
  }
# TO DO FIX THIS MESS

  if(verbose){
      print(paste0(bin_path," callVariants --refFile=",ref_genome, paste0(" --bamFiles=",tumor,",",norm)," --source=",vcf_overlay," --output=",out_file," --filterReadPairsWithSmallInserts=0 --minPosterior=0 --getVariantsFromBAMs=1  --logFileName=",paste0(out_file,".log"),opt," --nCPU=",threads))

  }
  system(paste0(bin_path," callVariants --refFile=",ref_genome, paste0(" --bamFiles=",tumor,",",norm), " --source=",vcf_overlay," --output=",out_file," --filterReadPairsWithSmallInserts=0 --minPosterior=0 --getVariantsFromBAMs=1 --logFileName=",paste0(out_file,".log"),opt," --nCPU=",threads))
  system(paste("cp", out_file, paste0(out_file,".tmp")))
  bgzip(bin_path=bin_path2,file=out_file)
  tab_indx(bin_path=bin_path3,file=paste0(out_file,".gz"))
  system(paste("cp", paste0(out_file,".tmp"), out_file))
  system(paste("rm -rf", paste0(out_file,".tmp")))
}
