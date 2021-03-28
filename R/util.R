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
#' @param vcf [Required] Path to VCF file to split.
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

  name=system(paste(bin_path," query -l ",vcf), intern = TRUE, ignore.stderr = TRUE)
  sapply(name,function(x){system(paste(bin_path, "view -c1 -Oz -s",x, "-o", paste0(out_file,"/",x,".vcf.gz"),vcf))
  })

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

  if (qual!=""){
    qual=paste0("\'%QUAL>",qual)
  }else{
    if (filter!=""){
      filter=paste0(" FILTER=\"",filter,"\" ")
    }
  }
  if (filter!=""){
    filter=paste0(" & FILTER=\"",filter,"\" ")
  }

  if (mq!=""){
    mq=paste0("& MQ>",mq,"\'")
  }



  if(verbose){
    print(paste(bin_path,"view  -i ",qual,filter,state,type,ref,mq,unfil_vcf,">",out_file))
  }
  system(paste(bin_path,"view  -i ",qual,filter,state,type,ref,mq,unfil_vcf,">",out_file))
  system(paste("cp", out_file, paste0(out_file,".tmp")))
  bgzip(bin_path=bin_path2,file=out_file)
  tab_indx(bin_path=bin_path3,file=paste0(out_file,".gz"))
  system(paste("cp", paste0(out_file,".tmp"), out_file))
  system(paste("rm -rf", paste0(out_file,".tmp")))
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



#' VCF stats file concatenation
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


#' This function filter and formats heterozygous SNP data for downstream analysis using clonet
#'
#' This function takes a path of a directory of unfiltered VCFs files generated using Platypus and
#' the a path of the directory of BAMs from which those VCFs were generated.
#'
#' @param bin_path Path to bcftools binary. Default tools/bcftools/bcftools.
#' @param bin_path2 Path to bgzip binary. Default tools/htslib/bgzip.
#' @param bin_path3 Path to tabix binary. Default tools/htslib/tabix.
#' @param bin_path4 Path to ASEQ binary. Default tools/ASEQ/binaries/linux64/ASEQ
#' @param unfil_vcf_dir Path to unfiltered VCF file directory.
#' @param qual Quality filter. Default 30.
#' @param mq Mapping quality filter. Default 40.
#' @param min_cov Minimum coverage to filter. Default 20.
#' @param bam_dir Path to BAM file directory.
#' @param patient_id Patient identifier.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads to use. Default 3
#' @export

format_SNP_data=function(bin_path="tools/bcftools/bcftools",bin_path2="tools/htslib/bgzip",bin_path3="tools/htslib/tabix",bin_path4="tools/ASEQ/binaries/linux64/ASEQ",unfil_vcf_dir="",bam_dir="",qual=30,mq=40,min_cov=20,patient_id="",verbose=FALSE,output_dir="",threads=3){
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
    files2=list.files(out_file_dir,recursive=TRUE,full.names=TRUE,pattern="FILTERED")
    files2=files2[grepl("vcf$",files2)]
    files2=as.data.frame(files2)
    names(files2)="VCF_path"
    files2$Sample=apply(files2,1,FUN=ULPwgs::get_sample_name)

    files3=list.files(bam_dir,recursive=TRUE,full.names=TRUE,pattern=patient_id)
    files3=files3[grepl("bam$",files3)]
    files3=as.data.frame(files3)
    names(files3)="BAM_path"
    files3$Sample=apply(files3,1,FUN=ULPwgs::get_sample_name)
    files=dplyr::left_join(files2,files3,by="Sample")
    print(files)
    pbapply::pblapply(X=1:nrow(files), FUN=function(x){call_ASEQ(vcf=as.character(files[x,1]),bin_path=bin_path4,bam=as.character(files[x,3]),mrq=mq,mbq=qual,mdc=min_cov,output_dir=out_file_dir,threads=1,verbose=verbose)},cl=cl)
    files3=list.files(out_file_dir,recursive=TRUE,full.names=TRUE,pattern="PILEUP.ASEQ")
    pbapply::pbapply(X=as.data.frame(files3),1,FUN=format_ASEQ_pileup,verbose=verbose,output_dir=out_file_dir,cl=cl)
    on.exit(parallel::stopCluster(cl))
  }
