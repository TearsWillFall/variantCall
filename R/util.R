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


#' Get SV type from svaba generated VCF
#'
#' As given by sfrenk comment on https://github.com/walaj/svaba/issues/4 annotates SV type
#' in svaba vcf file
#'
#' @param vcf svaba generated vcf to annotate. Get vcf
#' @export

annotate_sv_type <- function(vcf=""){
    # Find mate pair
    cols <- system(paste0('grep -v "##" ', vcf,' | grep "#" | sed s/#//'),intern=TRUE)
    cols <- strsplit(cols,"\t")[[1]]
    svaba_uniq <<- read.table(vcf, col.names = cols, stringsAsFactors = FALSE)
    svaba_uniq$SVTYPE <<- sapply(svaba_uniq$ID, FUN=get_sv_type,dat=svaba_uniq)
    fil=file(paste0(ULPwgs::get_sample_name(vcf),".svaba.sv.annotated.vcf"))
    writeLines(system(paste0(' grep "##" ', vcf ),intern=TRUE),fil)
    writeLines(paste0("#",paste0(cols,collapse="\t")),fil)
    write.table(svaba_uniq,append=TRUE,quote=FALSE,col.names=FALSE)
    close(fil)
}

#' Get SV type from svaba generated VCF
#'
#' As given by sfrenk comment on https://github.com/walaj/svaba/issues/4 annotates SV based on its breakpoints
#'
#'
#' @param x svaba generated vcf to annotate. Get vcf
#' @param dat Table with VCF data
#' @export

get_sv_type <- function(x,dat){
  root <- gsub(":[12]", "", x)
  mate1 <- paste0(root, ":1")
  mate2 <- paste0(root, ":2")
  print(dat)
  alt1 <- dat %>% filter(ID == mate1) %>% .$ALT
  alt2 <- dat %>% filter(ID == mate2) %>% .$ALT
  # Determine sv type based on breakpoint orientation
  if ((grepl("\\[", alt1) & grepl("\\[", alt2)) | (grepl("\\]", alt1) & grepl("\\]", alt2))){
      sv_type <- "INV"

  } else if (grepl("[A-Z]\\[", alt1) & grepl("^\\]", alt2)){
      sv_type <- "DEL"

  } else if (grepl("^\\]", alt1) & grepl("[A-Z]\\[", alt2)){
      sv_type <- "DUP/INS"

  } else{
      sv_type <- "UNKNOWN"
  }

  return(sv_type)
}


#' Learn Read Orientation Model (Mutect2)
#'
#' This function takes a f1r2 read information and generates an orientation model
#' to be used for further filtering steps.
#'
#' @param bin_path [REQUIRED] Path to GATK binary. Default tools/htslib/tabix.
#' @param f1r2 [OPTIONAL] Path to f1r2 file. Only if f1r2_dir is not given.
#' @param f1r2_dir [OPTIONAL] Path to f1r2 file dir. Only if f1r2 is not given.
#' @param verbose [Optional] Enables progress messages. Default False.
#' @export

learn_orientation=function(bin_path="tools/gatk/gatk",f1r2="",f1r2_dir="",output_name="",output_dir="",verbose=FALSE){

  if (f1r2!=""){
    f1r2=paste0(" -I ",paste0(f1r2,collapse=" -I "))
  }else{
    files=list.files(f1r2_dir,full.names=TRUE)
    files=files[grepl(".f1r2.tar.gz$",files)]
    f1r2=paste0(" -I ",paste0(files,collapse=" -I "))
  }

  if (output_name==""){
    sample_name=ULPwgs::get_sample_name(f1r2[1])
  }else{
    sample_name=output_name
  }

  sep="/"

  if(output_dir==""){
    sep=""
  }
  out_file_dir=paste0(output_dir,sep,sample_name,"_ORIENTATION_MODEL")

  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  if (verbose){

    print(paste0(bin_path," LearnReadOrientationModel ",f1r2, " -O ", paste0(out_file_dir,"/",sample_name,".read-orientation-model.tar.gz")))
  }
  system(paste0(bin_path," LearnReadOrientationModel ",f1r2, " -O ", paste0(out_file_dir,"/",sample_name,".read-orientation-model.tar.gz")))
}


#' Generate a set between a group of vcfs usign bcftools
#'
#' This function takes a path to either a vcf file or to a directory of vcf files
#' and outputs an inset/outset between these vcf files.
#'
#' @param bin_path [REQUIRED] Path to bcftools binary. Default tools/bcftools/bcftools
#' @param set_formula [REQUIRED] Formula to filter variants. Default =+2.
#' @param vcf [OPTIONAL] Path to a single vcf file or a vector of vcf files. Only if vcf_dir is not given.
#' @param vcf_dir [OPTIONAL] Path to a directory with vcf files to generate a set for. Only if vcf is not given.
#' @param filter [OPTIONAL] Filter variants by. Default PASS
#' @param output_dir [OPTIONAL] Path to output directory. Default current directory
#' @param verbose [Optional] Enables progress messages. Default False.
#' @export

vcf_sets=function(bin_path="tools/bcftools/bcftools",vcf="",vcf_dir="",set_formula="=+2",filter="PASS",output_dir="",verbose=FALSE){

  sep="/"

  if(output_dir==""){
    sep=""
    output_dir="SET"
  }

  if (filter!=""){
    filter=paste0(" -f ",filter)
  }

  if (!dir.exists(output_dir)){
      dir.create(output_dir,recursive=TRUE)
  }

  if (vcf_dir!=""){
    vcfs=list.files(vcfs,full.names=TRUE)
    vcfs=vcfs[grepl(".vcf.gz$",files)]
    vcfs=paste0(vcfs,collapse=" ")
  }else{
    vcfs=paste0(" ",vcf,collapse=" ")
  }

  private=paste0(" -n",set_formula)

  if (verbose){

    print(paste0(bin_path," isec -p ", output_dir,vcfs,filter,private))

  }
  system(paste0(bin_path," isec -p ", output_dir,vcfs,filter,private))
}



#' Generate a 1:n sets between a group of vcfs
#'
#' This function takes a path to either a vcf file or to a directory of vcf files
#' and outputs an inset/outset between these vcf files.
#'
#' @param bin_path [REQUIRED] Path to bcftools binary. Default tools/bcftools/bcftools
#' @param vcf [OPTIONAL] Path to a single vcf file or a vector of vcf files. Only if vcf_dir is not given.
#' @param vcf_dir [OPTIONAL] Path to a directory with vcf files to generate a set for. Only if vcf is not given.
#' @param filter [OPTIONAL] Filter variants by. Default PASS
#' @param output_dir [OPTIONAL] Path to output directory. Default current directory
#' @param set_names [OPTIONAL] Names to use to identify each set.
#' @param plot [OPTIONAL] Generate a VennDiagram plot with set intersections
#' @param threads [OPTIONAL] Number of threads to use. Default 3
#' @param verbose [Optional] Enables progress messages. Default False.
#' @export

generate_sets=function(bin_path="tools/bcftools/bcftools",vcf="",vcf_dir="",filter="PASS",output_dir="",set_names="",verbose=FALSE,threads=3,plot=TRUE){

  sep="/"
  if(output_dir==""){
    sep=""
  }

  out_file_dir=paste0(output_dir,sep,"SETS")

  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir,recursive=TRUE)
  }


  if (vcf_dir!=""){
    vcfs=list.files(vcfs,full.names=TRUE)
    vcfs=vcfs[grepl(".vcf.gz$",files)]
    n_vcfs=length(vcfs)
    vcfs=paste0(vcfs,collapse=" ")
  }else{
    n_vcfs=length(vcf)
    vcfs=paste0(" ",vcf,collapse=" ")
  }

  parallel::mclapply(X=1:n_vcfs,FUN=function(x){
    if (x>1 & x<n_vcfs){
        n <- n_vcfs
        m <- expand.grid(rep(list(0:1),n))
        m <- m[rowSums(m)==x ,]
        lapply(X=1:nrow(m),FUN=function(y){
          vcf_sets(bin_path=bin_path,vcf=vcf,vcf_dir=vcf_dir,set_formula=paste0("~",paste0(m[y,],collapse="")),filter=filter,output_dir=paste0(out_file_dir,"/SET_",x,"/SET_",paste0(m[y,],collapse="")),verbose=verbose)
        })
    }else{
      vcf_sets(bin_path=bin_path,vcf=vcf,vcf_dir=vcf_dir,set_formula=paste0("=",x),filter=filter,output_dir=paste0(out_file_dir,"/SET_",x),verbose=verbose)
    }
  })
  files=list.files(out_file_dir,recursive=TRUE,full.names=TRUE)
  files=files[grepl("sites.txt",files)]
  tables=lapply(files,FUN=read.table,colClasses="character")
  tables=dplyr::bind_rows(tables)
  tables=tables %>% dplyr::mutate(Variant=paste(V1,V2,V3,V4,sep="_"))
  colnames(tables)=c("chr","pos","ref","alt","inter","ID")
  sets=lapply(X=1:n_vcfs,FUN=function(z){
    search=rep(".",n_vcfs);
    search[z]="1";
    return(tables[grepl(paste0(search,collapse=""),tables$inter),]$ID);
  })

  if(set_names==""){
    names(sets)=paste0("Set_",1:n_vcfs)
  }else{
    names(sets)=set_names
  }
  if (plot){
    p1=ggvenn::ggvenn(sets)
    p1=p1+ggplot2::labs(title="Variant sets",subtitle=paste0("N=",nrow(tables)))
    ggplot2::ggsave(paste0(out_file_dir,"/variantSets_VennDiagram.png"),height=10,width=10)
  }
  write.table(tables,file=paste0(out_file_dir,"/variantData.txt"),quote=FALSE,row.names=FALSE,col.names=TRUE)
}

#' Filter Variant Tranches (Mutect2)
#'
#' This function takes an annotated VCF with CNN_D1 or CNN_D2 scores and filter the
#' variants based on set threshold.
#'
#' @param bin_path [REQUIRED] Path to GATK binary. Default tools/gatk/gatk
#' @param vcf [REQUIRED] Path to annotated VCF with CNN_D1/CNN_D2 scores
#' @param resources [OPTIONAL] Path to resources for variant filtering
#' @param output_name [OPTIONAL] Name of the sample to output
#' @param info_key [OPTIONAL] Annotation column to select. Default CNN_D1
#' @param snp_tranche [OPTIONAL] SNP tranche filter value. Default 99.95
#' @param indel_tranche [OPTIONAL] Indel tranche filter value. Default 99.4
#' @param keep_previous_filters [OPTIONAL] Keep previous filters in VCF. Default False
#' @param output_dir [OPTIONAL] Path to output dir
#' @param verbose [Optional] Enables progress messages. Default False
#' @export

FilterVariantTranches=function(bin_path="tools/gatk/gatk",vcf="",resources="",output_name="",info_key="CNN_1D",snp_tranche=99.95,indel_tranche=99.4,output_dir="",keep_previous_filters=FALSE,verbose=FALSE){

  prev_filters=" --invalidate-previous-filters "
  if (keep_previous_filters){
    prev_filters=" "
  }

  if (output_name==""){
    sample_name=ULPwgs::get_sample_name(vcf[1])
  }else{
    sample_name=output_name
  }

  if(resources!=""){
    resources=paste0(" --resource ",paste0(resources,collapse=" --resource "))
  }

  sep="/"

  if(output_dir==""){
    sep=""
  }
  out_file_dir=paste0(output_dir,sep,sample_name,"_FILTERED_TRENCHES")

  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }


  if (verbose){
    print(paste0(bin_path," FilterVariantTranches -V ",vcf," --info-key ",info_key, " -O ", paste0(out_file_dir,"/",sample_name,".FILTERED.vcf.gz")," --snp-tranche ",snp_tranche," --indel-tranche ",indel_tranche,resources,prev_filters))
  }
  system(paste0(bin_path," FilterVariantTranches -V ",vcf, " --info-key ",info_key, " -O ", paste0(out_file_dir,"/",sample_name,".FILTERED.vcf.gz")," --snp-tranche ",snp_tranche," --indel-tranche ",indel_tranche,resources,prev_filters))
}



#' Annotate a VCF with scores from a Convolutional Neural Network (CNN) {Mutect2}
#'
#' This function takes a vcf and then annotates it with CNN scores.
#' If a BAM files is supplied a 2D model will be applied, otherwise a 1D model will be applied.
#'
#' @param bin_path [REQUIRED] Path to GATK binary. Default tools/gatk/gatk
#' @param vcf [REQUIRED] Path to VCF file to annotate.
#' @param ref_genome [REQUIRED] Path to reference genome.
#' @param bam [OPTIONAL] Path to BAM file.
#' @param output_dir [OPTIONAL] Path to BAM file.
#' @param output_name [OPTIONAL] Default output name. If not supplied it will use the vcf file name.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @export


CNNScoreVariants=function(bin_path="tools/gatk/gatk",vcf="",ref_genome="",bam="",output_dir="",output_name="",verbose=FALSE){

  sep="/"
  if(output_dir==""){
    sep=""
  }

  sample_name=ULPwgs::get_sample_name(vcf)

  if(output_name!=""){
    sample_name=output_name
  }

  out_file_dir=paste0(output_dir,sep,sample_name,"_CNNscored")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  opt=""
  if (bam!=""){
    bam=paste0(" -I ",bam)
    opt=" -tensor-type read_tensor "
    out_file=paste0(out_file_dir,"/",sample_name,".CNNscored.2D.vcf")
  }else{
    out_file=paste0(out_file_dir,"/",sample_name,".CNNscored.1D.vcf")
  }

  if(verbose){
    print(paste0("source activate gatk; ",bin_path," CNNScoreVariants -V ",vcf," -R ",ref_genome, " -O ", out_file,bam,opt))
  }
  system(paste0("source activate gatk; ",bin_path," CNNScoreVariants -V ",vcf," -R ",ref_genome, " -O ", out_file,bam,opt))

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
        dir.create(out_file,recursive=TRUE)
    }

  if (verbose){
    name=system(paste(bin_path," query -l ",vcf), intern = TRUE, ignore.stderr = TRUE)
    print(paste(bin_path," query -l ",vcf))
    sapply(name,function(x){print(paste(bin_path, "view -c1 -Oz -s",x, "-o",paste0(out_file,"/",x,".vcf.gz"),vcf));
      system(paste(bin_path, "index", paste0(out_file,"/",x,".vcf.gz"),vcf))
    })
  }
  name=system(paste(bin_path," query -l ",vcf), intern = TRUE, ignore.stderr = TRUE)
  sapply(name,function(x){system(paste(bin_path, "view -c1 -Oz -s",x, "-o", paste0(out_file,"/",x,".vcf.gz"),vcf));
        system(paste(bin_path, "index", paste0(out_file,"/",x,".vcf.gz"),vcf))
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
      dir.create(out_file_dir,recursive=TRUE)
  }

  out_file=paste0(out_file_dir,"/",sample_name,".FILTERED.vcf")

  if (qual!=""){
    qual=paste0("\'%QUAL>",qual)
  }else{
    if (mq!=""){
      mq=paste0("\'%MQ>",mq)
    }else{
        if (filter!=""){
          filter=paste0("\'%FILTER=\"",filter,"\" ")
        }else{
          if (type!=""){
            type=paste0("\'%TYPE=\"",type,"\" ")
          }else{
            if (ref!=""){
              ref=paste0("\'%ID!=\"",ref,"\" ")
            }else{
              if (state!=""){
                state=paste0("\'%GT[0]=\"",state,"\" ")
              }
            }
          }
        }
      }
    }


    if (state!=""&!grepl("%",state)){
      state=paste0(" & GT[0]=\"",state,"\" ")
    }

    if (ref!=""&!grepl("%",ref)){
      ref=paste0(" & ID!=\"",ref,"\" ")
    }

    if (type!=""&!grepl("%",type)){
      type=paste0(" & TYPE=\"",type,"\" ")
    }
    if (filter!=""&!grepl("%",filter)){
      filter=paste0(" & FILTER=\"",filter,"\" ")
    }
    if (mq!=""&!grepl("%",mq)){
      mq=paste0(" & MQ>",mq)
    }

  if(verbose){
    print(paste(bin_path,"view  -i ",qual,filter,state,type,ref,mq,"\'",unfil_vcf,">",out_file))
  }
  system(paste(bin_path,"view  -i ",qual,filter,state,type,ref,mq,"\'",unfil_vcf,">",out_file))
  system(paste("cp", out_file, paste0(out_file,".tmp")))
  bgzip(bin_path=bin_path2,file=out_file)
  tab_indx(bin_path=bin_path3,file=paste0(out_file,".gz"))
  system(paste("cp", paste0(out_file,".tmp"), out_file))
  system(paste("rm -rf", paste0(out_file,".tmp")))
}


#' VCF file concatenation
#'
#' This function concatenates VCF files found in a directory.
#'
#' @param bin_path Path to bcftools binary. Default tools/bcftools/bcftools.
#' @param vcf Path to vcf files to concatenate. Only if vcf_dir not given
#' @param vcf_dir Path to directory with vcf files to concatenate. Only if vcf not given
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export



vcf_concatenate=function(bin_path="tools/bcftools/bcftools",vcf="",vcf_dir="",verbose=FALSE,output_dir=""){
  sep="/"
  if(output_dir==""){
    sep=""
  }
  if (vcf_dir!=""){
    files=list.files(vcf_dir,full.names=TRUE)
    files=files[grepl(".vcf$",files)]
    sample_name=ULPwgs::get_sample_name(list.files(vcf_dir)[1])
  }else{
    files=vcf
    sample_name=ULPwgs::get_sample_name(vcf[1])
  }

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

#' Filter VEP results
#'
#' This function filters VCFs generated by VEP
#'
#' @param bin_path Path to vep binary. Default tools/ensembl-vep/filter_vep.
#' @param bin_path2 Path to bgzip binary. Default tools/htslib/bgzip.
#' @param bin_path3 Path to tabix binary. Default tools/htslib/tabix.
#' @param unf_vcf Path to vcf file.
#' @param filter Column tags to filter by. Examples "(MAX_AF < 0.01 or not MAX_AF) and FILTER = PASS and SYMBOL in /home/regmova/PCFselect/Panel_V2/PANEL_GENE.txt". For more information check VEP filter info.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export

filter_VEP=function(bin_path="tools/ensembl-vep/filter_vep",bin_path2="tools/htslib/bgzip",bin_path3="tools/htslib/tabix",unf_vcf="",filter="",verbose=FALSE,output_dir=""){
  sep="/"
  if(output_dir==""){
    sep=""
  }
  sample_name=ULPwgs::get_sample_name(unf_vcf)
  out_file_dir=paste0(output_dir,sep,sample_name,"_FILTERED_VEP")
  if (!dir.exists(out_file_dir)){
      dir.create(out_file_dir)
  }

  out_file=paste0(out_file_dir,"/",sample_name,".FILTERED.VEP.vcf")

  if(verbose){
    print(paste(bin_path,"-i ",unf_vcf,"-f ",filter,">",out_file))
  }
  system(paste(bin_path,"-i ",unf_vcf,"-f ",filter,">",out_file))
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
#' @param worst [OPTIONAL]  Only keep the variants with the worst outcome.
#' @param output_name [OPTIONAL]  Name of output file.
#' @param verbose [OPTIONAL]  Enables progress messages. Default False.
#' @export


format_PM_data=function(bin_path="tools/bcftools/bcftools",vcf_dir="",pattern="",vcf="",worst=FALSE,output_dir="",output_name="",verbose=FALSE){

  if(vcf!="" & vcf_dir!=""){

    stop("Only vcf or vcf_dir arguments can be provided, not both.")
  }

  if (vcf_dir!=""){
    files=list.files(path=vcf_dir,pattern=pattern,full.names=TRUE,recursive=TRUE)
    files=files[grepl("vcf$",files)]
  }else {
    files=vcf
  }

  if(output_name==""){
    output_name="PMreadsCount"
  }
  wrst=""
  if(worst){
    wrst= " -s worst "
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
    lapply(files,FUN=function(x){print(paste0(bin_path," +", system(paste0("echo ",sub("bcftools$", "plugins\\1",bin_path),"/split-vep.so"),intern=TRUE)," -f \'%SYMBOL{0}\\t%CHROM\\t%POS\\t%FILTER\\t[\\t%AD]\\n\' -i \'TYPE=\"snp\" & AD[0:1]>0\' ",wrst, x," | awk -F \'[\t,]\' \'{ print  $1,$2,$3,$3, \"",ULPwgs::get_sample_name(x),"\", $5, $6 }\' OFS=\'\\t\' >> ",out_file))})
    print(paste("printf","\'Gene.id\tchr\tpos\tref\talt\tFILTER\tAD\tDP\tAllele\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature_type\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tDISTANCE\tSTRAND\tFLAGS\tVARIANT_CLASS\tSYMBOL_SOURCE\tHGNC_ID\tCANONICAL\tMANE\tTSL\tAPPRIS\tCCDS\tENSP\tSWISSPROT\tTREMBL\tUNIPARC\tUNIPROT_ISOFORM\tGENE_PHENO\tSIFT\tPolyPhen\tDOMAINS\tmiRNA\tHGVS_OFFSET\tAF\tAFR_AF\tAMR_AF\tEAS_AF\tEUR_AF\tSAS_AF\tAA_AF\tEA_AF\tgnomAD_AF\tgnomAD_AFR_AF\tgnomAD_AMR_AF\tgnomAD_ASJ_AF\tgnomAD_EAS_AF\tgnomAD_FIN_AF\tgnomAD_NFE_AF\tgnomAD_OTH_AF\tgnomAD_SAS_AF\tMAX_AF\tMAX_AF_POPS\tCLIN_SIG\tSOMATIC\tPHENO\tPUBMED\tMOTIF_NAME\tMOTIF_POS\tHIGH_INF_POS\tMOTIF_SCORE_CHANGE\tTRANSCRIPTION_FACTORS\tsample\n\'",">",paste0(out_file,".ALL")))
    lapply(files,FUN=function(x){print(paste0(bin_path," +",system(paste0("echo ",sub("bcftools$", "plugins\\1",bin_path),"/split-vep.so"),intern=TRUE), " -f \'%SYMBOL{0}\\t%CHROM\\t%POS\\t%FILTER\\t[\\t%AD]\\t[\\t%DP]\\t%CSQ\\n\' -i \'AD[0:1]>0\' ",wrst," -d -A tab ",x," |sed \'s/$/\\t",ULPwgs::get_sample_name(x),"/\' >>",paste0(out_file,".ALL")))})
  }
    system(paste("printf","\'Gene.id\tchr\tstart\tend\tsample\tref.count\talt.count\n\'",">",out_file))
    lapply(files,FUN=function(x){system(paste0(bin_path," +",system(paste0("echo ",sub("bcftools$", "plugins\\1",bin_path),"/split-vep.so"),intern=TRUE), " -f \'%SYMBOL{0}\\t%CHROM\\t%POS\\t%FILTER\\t[\\t%AD]\\n\' -i \'TYPE=\"snp\" & AD[0:1]>0\' ",wrst ,x," | awk -F \'[\t,]\' \'{ print  $1,$2,$3,$3, \"",ULPwgs::get_sample_name(x),"\", $5, $6 }\' OFS=\'\\t\' >> ",out_file))})

    system(paste("printf","\'Gene.id\tchr\tpos\tref\talt\tFILTER\tAD\tDP\tAllele\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature_type\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tDISTANCE\tSTRAND\tFLAGS\tVARIANT_CLASS\tSYMBOL_SOURCE\tHGNC_ID\tCANONICAL\tMANE\tTSL\tAPPRIS\tCCDS\tENSP\tSWISSPROT\tTREMBL\tUNIPARC\tUNIPROT_ISOFORM\tGENE_PHENO\tSIFT\tPolyPhen\tDOMAINS\tmiRNA\tHGVS_OFFSET\tAF\tAFR_AF\tAMR_AF\tEAS_AF\tEUR_AF\tSAS_AF\tAA_A\tEA_AF\tgnomAD_AF\tgnomAD_AFR_AF\tgnomAD_AMR_AF\tgnomAD_ASJ_AF\tgnomAD_EAS_AF\tgnomAD_FIN_AF\tgnomAD_NFE_AF\tgnomAD_OTH_AF\tgnomAD_SAS_AF\tMAX_AF\tMAX_AF_POPS\tCLIN_SIG\tSOMATIC\tPHENO\tPUBMED\tMOTIF_NAME\tMOTIF_POS\tHIGH_INF_POS\tMOTIF_SCORE_CHANGE\tTRANSCRIPTION_FACTORS\tsample\n\'",">",paste0(out_file,".ALL")))
    lapply(files,FUN=function(x){system(paste0(bin_path," +",system(paste0("echo ",sub("bcftools$", "plugins\\1",bin_path),"/split-vep.so"),intern=TRUE), " -f \'%SYMBOL{0}\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t[\\t%AD]\\t[\\t%DP]\\t%CSQ\\n\' -i \' AD[0:1]>0\' ",wrst," -d -A tab ",x," |sed \'s/$/\\t",ULPwgs::get_sample_name(x),"/\' >>",paste0(out_file,".ALL")))})
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
#' @param bin_path [REQUIRED] Path to gatk binary. Default tools/gatk/gatk.
#' @param bin_path2 [REQUIRED] Path to bgzip binary. Default tools/htslib/bgzip.
#' @param bin_path3 [REQUIRED] Path to tabix binary. Default tools/htslib/tabix.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param unfil_vcf [REQUIRED] Path to unfiltered vcf file.
#' @param unfil_vcf_stats [REQUIRED] Path to unfiltered vcf file stats.
#' @param contamination [OPTIONAL] Path to contamination table. Also requires segmentation.
#' @param segmentation [OPTIONAL] Path to segments table. Also requires contamination.
#' @param orientation_model [OPTIONAL] Path to orientation model generated using learn_orientation.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @export


vcf_filtering=function(bin_path="tools/gatk/gatk",bin_path2="tools/htslib/bgzip",bin_path3="tools/htslib/tabix",unfil_vcf="",ref_genome="",unfil_vcf_stats="",contamination="",segmentation="",orientation_model="",verbose=FALSE,output_dir=""){

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

    if (contamination!="" & segmentation !=""){
      contamination=paste0(" --contamination-table ", paste0(contamination,collapse=" --contamination-table "))
      segmentation=paste0(" --tumor-segmentation ",paste0(segmentation,collapse=" --tumor-segmentation "))
    }

    if (orientation_model!=""){
      orientation_model=paste0(" --ob-priors ",orientation_model)
    }

    if(verbose){
      print(paste(bin_path,"FilterMutectCalls -O",out_file," -R ",ref_genome," -V ",unfil_vcf," -stats ",unfil_vcf_stats, contamination,segmentation,orientation_model))
    }
    system(paste(bin_path,"FilterMutectCalls -O",out_file," -R ",ref_genome," -V ",unfil_vcf," -stats ",unfil_vcf_stats, contamination,segmentation,orientation_model))
    system(paste("cp", out_file, paste0(out_file,".tmp")))
    bgzip(bin_path=bin_path2,file=out_file)
    tab_indx(bin_path=bin_path3,file=paste0(out_file,".gz"))
    system(paste("cp", paste0(out_file,".tmp"), out_file))
    system(paste("rm -rf", paste0(out_file,".tmp")))
  }

#' VCF contamination estimation using GATK
#'
#' This function estimates sample contamination for further variant filtering
#'
#' @param bin_path [REQUIRED] Path to gatk binary. Default tools/gatk/gatk.
#' @param pileup_table_tumor [REQUIRED] Path to pileup_table for tumor sample generated using pileup_summary
#' @param pileup_table_normal [OPTIONAL] Path to pileup_table for matched normal sample generated using pileup_summary
#' @param output_dir [OPTIONAL]  Path to the output directory.
#' @param verbose [OPTIONAL]  Enables progress messages. Default False.
#' @export


estimate_contamination=function(pileup_table_tumor="",bin_path="tools/gatk/gatk",pileup_table_normal="",verbose=FALSE,output_dir=""){

      sep="/"
      if(output_dir==""){
        sep=""
      }
      sample_name=ULPwgs::get_sample_name(pileup_table_tumor)
      if (!dir.exists(output_dir)){
          dir.create(output_dir)
      }

      out_file=paste0(output_dir,"/",sample_name,".contamination.table")
      out_file2=paste0(output_dir,"/",sample_name,".segments.table")
      if (pileup_table_normal!=""){
        pileup_table_normal=paste0("-matched ",pileup_table_normal)
      }

      if(verbose){
        print(paste(bin_path,"CalculateContamination -O ",out_file," -tumor-segmentation ",out_file2," -I ",pileup_table_tumor,pileup_table_normal))
      }
      system(paste(bin_path,"CalculateContamination -O ",out_file," -tumor-segmentation ",out_file2," -I ",pileup_table_tumor,pileup_table_normal))
  }

#' VCF contamination estimation using GATK
#'
#' This function estimates sample contamination for further variant filtering
#'
#' @param bin_path [REQUIRED] Path to gatk binary. Default tools/gatk/gatk.
#' @param bam_dir [REQUIRED] Path to directory with BAM files
#' @param germ_pattern [REQUIRED] Pattern to match normal samples.
#' @param patient_id [REQUIRED] Pattern to match patient specific samples.
#' @param db [REQUIRED] Path to vcf with known variants.
#' @param interval [REQUIRED] Path to vcf with intervals to analyze.
#' @param output_dir [OPTIONAL]  Path to the output directory.
#' @param verbose [OPTIONAL]  Enables progress messages. Default False.
#' @param threads [OPTIONAL]  Number of threads to use. Default 3.
#' @export


estimate_contamination_parallel=function(bin_path="tools/gatk/gatk",bam_dir="",germ_pattern="GL",patient_id="",db="",interval="",verbose=FALSE,output_dir="",threads=3){

      sep="/"
      if(output_dir==""){
        sep=""
      }
      sample_name=ULPwgs::get_sample_name(patient_id)
      out_file_dir=paste0(output_dir,sep,sample_name,"_CONTAMINATION")
      if (!dir.exists(out_file_dir)){
          dir.create(out_file_dir)
      }

      files=list.files(bam_dir,recursive=TRUE,full.names=TRUE,pattern=patient_id)
      files=files[grepl("bam$",files)]
      tumor_bams=files[!grepl(germ_pattern,files)]
      normal_bam=files[grepl(germ_pattern,files)]

      cl=parallel::makeCluster(threads)
      pbapply::pblapply(X=files,FUN=get_pileup_summary,bin_path=bin_path,db=db,interval=interval,verbose=verbose,output_dir=out_file_dir,cl=cl)
      on.exit(parallel::stopCluster(cl))

      files=list.files(out_file_dir,recursive=TRUE,full.names=TRUE,pattern=patient_id)
      files=files[grepl(".pileup.table$",files)]
      tumor_pileups=files[!grepl(germ_pattern,files)]
      normal_pileup=files[grepl(germ_pattern,files)]

      cl=parallel::makeCluster(threads)
      pbapply::pblapply(X=tumor_pileups,FUN=estimate_contamination,bin_path=bin_path,pileup_table_normal=normal_pileup,verbose=verbose,output_dir=out_file_dir,cl=cl)
      on.exit(parallel::stopCluster(cl))
}




#' BAM pileup summary for known sites using GATK
#'
#' This function generates a pileup summary that is further used in estimating contamination
#'
#' @param bin_path [REQUIRED] Path to gatk binary. Default tools/gatk/gatk.
#' @param bam [REQUIRED] Path to BAM file.
#' @param db [REQUIRED] Path to vcf with known variants.
#' @param interval [REQUIRED] Path to vcf with intervals to analyze.
#' @param output_dir [OPTIONAL]  Path to the output directory.
#' @param verbose [OPTIONAL]  Enables progress messages. Default False.
#' @export

get_pileup_summary=function(bam="",bin_path="tools/gatk/gatk",db="",interval="",verbose=FALSE,output_dir=""){
      sep="/"
      if(output_dir==""){
        sep=""
      }
      sample_name=ULPwgs::get_sample_name(bam)
      if (!dir.exists(output_dir)){
          dir.create(output_dir)
      }

      out_file=paste0(output_dir,"/",sample_name,".pileup.table")

      if(verbose){
        print(paste(bin_path,"GetPileupSummaries  -O",out_file," -V ",db," -L ",interval, " -I ",bam))
      }
      system(paste(bin_path,"GetPileupSummaries  -O",out_file," -V ",db," -L ",interval, " -I ",bam))

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

    pbapply::pblapply(X=1:nrow(files), FUN=function(x){call_ASEQ(vcf=as.character(files[x,1]),bin_path=bin_path4,bam=as.character(files[x,3]),mrq=mq,mbq=qual,mdc=min_cov,output_dir=out_file_dir,threads=1,verbose=verbose)},cl=cl)
    files3=list.files(out_file_dir,recursive=TRUE,full.names=TRUE,pattern="PILEUP.ASEQ")
    pbapply::pbapply(X=as.data.frame(files3),1,FUN=format_ASEQ_pileup,verbose=verbose,output_dir=out_file_dir,cl=cl)
    on.exit(parallel::stopCluster(cl))
  }
