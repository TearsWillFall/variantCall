#' Index tab separated genomic regions
#'
#' This function takes a .vcf.bgz tab separated genomic region file and generates an index for it
#'
#' @param bin_path [REQUIRED] Path to bcftools binary. Default tools/htslib/tabix.
#' @param file [REQUIRED] Path to VCF file to index.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
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
#' @param output_dir Path to directory to output results
#' @export

annotate_sv_type <- function(vcf="",output_dir=""){
  sep="/"
  if(output_dir==""){
    sep=""
  }

  out_file_dir=paste0(output_dir,sep)
  if (!dir.exists(out_file_dir) & !out_file_dir==""){
      dir.create(out_file_dir)
  }
  # Find mate pair

  cols <- system(paste0('grep -v "##" ', vcf,' | grep "#" | sed s/#//'),intern=TRUE)
  cols <- strsplit(cols,"\t")[[1]]
  svaba_uniq=tryCatch(
    {svaba_uniq = read.table(vcf, col.names = cols, stringsAsFactors = FALSE);
     svaba_uniq$INFO = paste0(svaba_uniq$INFO,";SVANOT=",sapply(svaba_uniq$ID, FUN=get_sv_type,dat=svaba_uniq));
     svaba_uniq}
        ,error=function(e){
      svaba_uniq=data.frame(matrix(ncol = length(cols), nrow = 0));
      colnames(svaba_uniq)=cols;
      return(svaba_uniq)})
  fil=paste0(out_file_dir,paste0(ULPwgs::get_sample_name(vcf),".svaba.sv.annotated.vcf"))
  cat(system(paste0('grep "##" ', vcf ),intern=TRUE),file=fil,sep="\n")
  cat('##INFO=<ID=SVANOT,Number=1,Type=String,Description=\"Structural variant annotation\">',file=fil,sep="\n",append=TRUE)
  cat(paste0("#",paste0(cols,collapse="\t")),file=fil,sep="\n",append=TRUE)
  write.table(x=svaba_uniq,file=fil,append=TRUE,quote=FALSE,col.names=FALSE,sep="\t",row.names=FALSE)
}



#' Compress and index vcf file
#'
#' This function takes an uncompressed vcf file and uses bgzip and tabix to
#' compress and index the file while preserving the original file
#'
#' @param vcf [REQUIRED] Path to VCF to compress
#' @param output_dir [OPTIONAL] Path to output directory
#' @export

compress_and_index_vcf=function(bin_path="tools/htslib/bgzip",bin_path2="tools/htslib/tabix",vcf="",output_dir=""){
  sep="/"
  if(output_dir==""){
    sep=""
  }
  out_file_dir=paste0(output_dir,sep)

  system(paste("cp",paste0(vcf), paste0(vcf,".tmp")))
  bgzip(bin_path=bin_path,file=vcf)
  tab_indx(bin_path=bin_path2,file=paste0(vcf,".gz"))
  system(paste("mv",paste0(vcf,".tmp"), vcf))
  if (!out_file_dir==""){
    if (!dir.exists(out_file_dir) ){
        dir.create(out_file_dir,recursive=TRUE)
    }
    system(paste("mv",paste0(vcf,".gz*"), out_file_dir))
  }
}


#' Intersect vcf with bed file
#'
#' Generate a set of variants covered by the regions found in/out the bed file
#'
#' @param vcf [REQUIRED] Path to VCF file
#' @param bed [REQUIRED] Path to BED file
#' @param output_name [OPTIONAL] Name for output_files. IF not given it will use input vcf name
#' @param output_dir [OPTIONAL] Path to output directory. if not given outputs to current directory.
#' @export

vcf_intersect_bed <- function(vcf="",bed="",output_name="",output_dir=""){
  sep="/"
  if(output_dir==""){
    sep=""
  }

  out_file_dir=paste0(output_dir,sep)
  if (!dir.exists(out_file_dir) & !out_file_dir==""){
      dir.create(out_file_dir)
  }
  if(output_name!=""){
    sample_name=output_name
  }else{
    sample_name=ULPwgs::get_sample_name(vcf)
  }
  cols <- system(paste0('grep -v "##" ', vcf,' | grep "#" | sed s/#//'),intern=TRUE)
  cols <- strsplit(cols,"\t")[[1]]
  variants = read.table(vcf, col.names = cols, stringsAsFactors = FALSE)
  variants$CHROM=as.character(variants$CHROM)
  regions= read.table(bed, stringsAsFactors = FALSE)
  regions$V1=as.character(sub("chr","",regions$V1))
  variants_in=dplyr::semi_join(variants,regions,by=c("CHROM"="V1","POS"="V3"))
  variants_out=dplyr::anti_join(variants,regions,by=c("CHROM"="V1","POS"="V3"))
  file_out=paste0(out_file_dir,sample_name,".IN_BED.vcf")
  file_out2=paste0(out_file_dir,sample_name,".OUT_BED.vcf")
  cat(system(paste0('grep "##" ', vcf ),intern=TRUE),file=file_out,sep="\n")
  cat(system(paste0('grep "##" ', vcf ),intern=TRUE),file=file_out2,sep="\n")
  cat(paste0('##Filter In-section between vcf=',vcf," and bed=", bed),file=file_out,sep="\n",append=TRUE)
  cat(paste0('##Filter Out-section between vcf=',vcf," and bed=", bed),file=file_out2,sep="\n",append=TRUE)
  cat(paste0("#",paste0(cols,collapse="\t")),file=file_out,sep="\n",append=TRUE)
  cat(paste0("#",paste0(cols,collapse="\t")),file=file_out2,sep="\n",append=TRUE)
  write.table(x=variants_in,file=file_out,append=TRUE,quote=FALSE,col.names=FALSE,sep="\t",row.names=FALSE)
  write.table(x=variants_out,file=file_out2,append=TRUE,quote=FALSE,col.names=FALSE,sep="\t",row.names=FALSE)
}


#' VCF communality
#'
#' Estimate % of common SNPs between two VCF files and returns the communality between them
#'
#' @param vcf [REQUIRED] Path to VCF file
#' @param vcf2 [REQUIRED] Path to VCF file
#' @param output_name [OPTIONAL] Name for output_files. IF not given it will use input vcf name
#' @param verbose [OPTIONAL] Extra verbose. Default FALSE
#' @param output_dir [OPTIONAL] Path to output directory. if not given outputs to current directory.
#' @export

vcf_communality =function(bin_path="tools/bcftools/bcftools",vcf="",vcf2="",output_name="",output_dir="",verbose=FALSE){

  sep="/"
  if(output_dir==""){
    sep=""
  }

  out_file_dir=paste0(output_dir,sep)
  if (!dir.exists(out_file_dir) & !out_file_dir==""){
      dir.create(out_file_dir)
  }
    # if(output_name!=""){
    #   sample_name=output_name
    # }else{
    #   sample_name=ULPwgs::get_sample_name(vcf)
    # }
  generate_sets(bin_path=bin_path,vcf=c(vcf,vcf2),filter="",output_dir=paste0(out_file_dir,paste0(ULPwgs::get_sample_name(vcf),"_AND_",ULPwgs::get_sample_name(vcf2))),verbose=verbose,set_names=c(ULPwgs::get_sample_name(vcf),ULPwgs::get_sample_name(vcf2)))
  not_common=as.numeric(system(paste0("cat ",paste0(out_file_dir,paste0(ULPwgs::get_sample_name(vcf),"_AND_",ULPwgs::get_sample_name(vcf2))),"/SETS/SET_1/sites.txt | wc -l"),intern=TRUE))
  common=as.numeric(system(paste0("cat ",paste0(out_file_dir,paste0(ULPwgs::get_sample_name(vcf),"_AND_",ULPwgs::get_sample_name(vcf2))),"/SETS/SET_2/sites.txt | wc -l"),intern=TRUE))
  communality=common/(common+not_common)
  return(data.frame(vcf=vcf,communality=communality,common=common,not_common))
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
  alt1 <- dat %>% dplyr::filter(ID == mate1) %>% .$ALT
  alt2 <- dat %>% dplyr::filter(ID == mate2) %>% .$ALT
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
  tables=lapply(files,FUN=function(x){
    tryCatch(read.table(x,colClasses="character"),error=function(x) NULL)})
  tables=tables[!sapply(tables,is.null)]
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
                state=paste0("\'GT[0]=\"",state,"\" ")
              }
            }
          }
        }
      }
    }


    if (state!=""&!grepl("GT",state)){
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

#' VCF fix MNPs
#'
#' Exclude MNPs from vcf file
#'
#' @param bin_path Path to bcftools binary. Default tools/bcftools/bcftools.
#' @param bin_path2 Path to bgzip binary. Default tools/htslib/bgzip.
#' @param bin_path3 Path to tabix binary. Default tools/htslib/tabix.
#' @param vcf Path to unfiltered vcf file.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export


fix_mnps=function(vcf="",bin_path="tools/bcftools/bcftools",bin_path2="tools/htslib/bgzip",bin_path3="tools/htslib/tabix",verbose=FALSE,output_dir=""){

  sep="/"
  if(output_dir==""){
    sep=""
  }
  sample_name=ULPwgs::get_sample_name(vcf)
  out_file_dir=paste0(output_dir,sep)
  if (!dir.exists(out_file_dir) & !out_file_dir==""){
      dir.create(out_file_dir,recursive=TRUE)
  }

  out_file=paste0(out_file_dir,sample_name,".vcf")

  if(verbose){
    print(paste(bin_path,"view --exclude-types mnps ",vcf,">",out_file))
  }
  system(paste(bin_path,"view --exclude-types mnps ",vcf,">",out_file))
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
      dir.create(out_file_dir,recursive=TRUE)
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
#' @param cols_to_keep [DEFAULT==c("chr","start","end","log2","sample"))] Columns to keep from original segmentation bed files
#' @param output_dir Path to the output directory.
#' @param output_name Name of the file to output.
#' @param verbose Enables progress messages. Default False.
#' @export



format_segment_data=function(seg_file="",dir_segment="",pattern="",cols_to_keep=c("chr","start","end","log2","sample"),output_dir="",output_name="",verbose=FALSE){

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
  names(data)[c(1,2,3,5,ncol(data))]=c("chr","start","end","log2","sample")
  data=data[,c(cols_to_keep)]

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
#' @param unfil_vcf Path to unfiltered VCF file.
#' @param unfil_vcf_dir Path to unfiltered VCF file directory.
#' @param qual Quality filter. Default 30.
#' @param mq Mapping quality filter. Default 40.
#' @param min_cov Minimum coverage to filter. Default 20.
#' @param bam_dir Path to BAM file directory.
#' @param patient_id Patient identifier.
#' @param germ_pattern Germline identifier.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default FALSE.
#' @param remove_tmp Remove temporary files. Default TRUE
#' @param threads Number of threads to use. Default 3
#' @export

format_SNP_data=function(bin_path="tools/bcftools/bcftools",bin_path2="tools/htslib/bgzip",bin_path3="tools/htslib/tabix",bin_path4="tools/ASEQ/binaries/linux64/ASEQ",unfil_vcf="",unfil_vcf_dir="",bam_dir="",germ_pattern="GL",qual=30,mq=40,min_cov=20,patient_id="",verbose=FALSE,output_dir="",threads=3,remove_tmp=TRUE){
    sep="/"
    if(output_dir==""){
      sep=""
    }
    out_file_dir=paste0(output_dir,sep,patient_id,"_PROCESSED_SNPs")
    if (!dir.exists(out_file_dir)){
        dir.create(out_file_dir)
    }
    cl=parallel::makeCluster(threads)
    if (unfil_vcf_dir!=""){
      files0=list.files(unfil_vcf_dir,recursive=TRUE,full.names=TRUE,pattern=patient_id)
      files0=files0[grepl("vcf.gz$",files0)]

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

    }else{
      files3=list.files(bam_dir,recursive=TRUE,full.names=TRUE,pattern=patient_id)
      files3=files3[grepl("bam$",files3)]
      pbapply::pblapply(X=1:length(files3), FUN=function(x){call_ASEQ(vcf=unfil_vcf,bam=as.character(files3[x]),mrq=mq,bin_path=bin_path4,mbq=qual,mdc=min_cov,output_dir=out_file_dir,threads=1,verbose=verbose)},cl=cl)
    }

    files3=list.files(out_file_dir,recursive=TRUE,full.names=TRUE,pattern="PILEUP.ASEQ")
    pbapply::pbapply(X=as.data.frame(files3),1,FUN=format_ASEQ_pileup,verbose=verbose,output_dir=out_file_dir,cl=cl)
    on.exit(parallel::stopCluster(cl))
    files=list.files(out_file_dir,recursive=TRUE,full.names=TRUE,pattern=".snp")
    tumor_snps=files[!grepl(germ_pattern,files)]
    germ_snps=files[grepl(germ_pattern,files)]
    out_file_dir2=paste0(out_file_dir,"/RESULTS")
    out_file_dir2_ger=paste0(out_file_dir2,"/NormalPileup")
    out_file_dir2_tumor=paste0(out_file_dir2,"/TumorPileup")
    if (!dir.exists(out_file_dir2_ger)){
        dir.create(out_file_dir2_ger,recursive=TRUE)
    }
    if (!dir.exists(out_file_dir2_tumor)){
        dir.create(out_file_dir2_tumor,recursive=TRUE)
    }

    system(paste("cp -rf -t ",out_file_dir2_tumor,paste(tumor_snps,collapse=" ")))
    system(paste("cp -rf -t ",out_file_dir2_ger,paste(germ_snps,collapse=" ")))
    if (remove_tmp){
      system(paste0("ls -d ", paste0(out_file_dir,"/*"), " |egrep -v RESULTS | xargs rm -rf"))
    }
  }


#' This function generates a config files for CLONET pipeline
#'
#' This function takes multiple parameters and generates a config file for the CLONET tool
## Stages to perform:
#   1. analyse single sample and produe RData
#   2. create beta table aggragating single samples analysis
#   3. compute ploidy shift
#   4. compute global admixture
#   5. compute clonality table
#   6. compute allele specific copy number table
#'
#' @param patient_id Patient ID. Default Patient
#' @param clonet_dir Path to CLONET dir
#' @param snp_dir Path to formated SNP data dir.
#' @param segment_data Path to file with formated segment data.
#' @param sample_info Path to file with sample info data.
#' @param min_snp_cov Minimum tumor coverage for informative SNPs. Default 10.
#' @param min_nsnps Minimum number of SNPs per segment. Default 10.
#' @param min_seg_cov Minimum segment coverage. Default 20.
#' @param equal_betaThr Minimum value of beta above which the two alleles are present in the same number. Default 0.9
#' @param max_homo_dels Homozygous deletions threshold. Default 0.01
#' @param del_log_thr Parameters of a valid deletion used to compute Adm.global. Default c(-1,-0.25)
#' @param alpha_par  Percentage of used deletions to compute Adm.global varibility interval. Default 0.9
#' @param clonal_thr Clonality value threshold. Default 0.85
#' @param beta_thr Beta value threshold. Default 0.85
#' @param stages Analysis stages to run trough.Default c(1,2,3,4,5,6)
#' @param comp_ref_map_bias Compute reference mapping bias and to adjust beta estimation. Default FALSE
#' @param beta_decimals Number of beta value decimals to report. Default 3
#' @param beta_method Method for beta estimation. Default STM. Options STM/GB
#' @param adm_method Method for admixture estimation. Default 2D. Options 1D/2D
#' @param jobs Number of Samples to analyze in parallel. Default 1
#' @param threads Number of threads per job to use. Default 3
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default FALSE.
#' @export

generate_CLONET_config=function(patient_id="PATIENT",clonet_dir="tools/CLONET",snp_dir="",segment_data="",
    sample_info="",min_snp_cov=10,min_nsnps=10,min_seg_cov=20,equal_betaThr=0.9,max_homo_dels=0.01,
    del_log_thr=c(-1,-0.25),alpha_par=0.9,clonal_thr=0.85,beta_thr=0.85,
    stages=c(1,2,3,4,5,6),comp_ref_map_bias=FALSE,beta_decimals=3,ale_imb_thr=0.5,
    beta_method="STM",adm_method="2D",jobs=1,threads=3,output_dir=""){

    sep="/"
    if(output_dir==""){
      sep=""
    }
    out_file_dir=paste0(output_dir,sep,patient_id,"_CLONET_CONFIG")

    if (!dir.exists(out_file_dir)){
        dir.create(out_file_dir,recursive=TRUE)
    }

    config_data=paste("
    ## Output directory
    output_DIR <-", paste0("'",output_dir,"/OUTPUT/'"),"

    ## Path to CLONET functions
    path_to_CLONET_functions <-", paste0("'",clonet_dir,"/'"),"

    # information about sample names
    sampleInfoFile <-",paste0("'",sample_info,"'"),"

    # list of segment
    segmentListFile <-", paste0("'",segment_data,"'"),"

    # folder with informative SNPs
    pileup_dir <- ",paste0("'",snp_dir,"/'"),"

    # Suffix of the informative SNPs pileup
    PaPI_Suffix <- '.snps'

    ## Path to error table
    errorTable_file =",paste0("'",clonet_dir,'/Examples/BasicData/errorTable.csv',"'"),"

    ## Method to compete beta
    # 'GB' : classic model defined in http://www.genomebiology.com/2014/15/8/439 (only for backward compatibility)
    # 'STM' : new method based on error estimation on matched normal published in http://stm.sciencemag.org/content/6/254/254ra125.short (RECOMMENDED)
    betaCompute.method <-",paste0("'",beta_method,"'"),"

    ## Method to compete Adm.global
    # '1D' : classic model that use only deletions (only for backward compatibility)
    # '2D' : model based on cnA vs cnB transformation (RECOMMENDED)
    adm.global.method <- ",paste0("'",adm_method,"'"),"

    # minimum tumor coverage to consider an informative SNP as valid
    minCoverage <-",min_snp_cov,"

    # number of samples to process in parallel
    NsamplesToProcessInParallel <-", jobs,"

    # number of cores assigned to the analysis of each process
    perSampleCores <-", threads,"

    # min number of informative SNPs for a genomic segment
    min_nsnps <-",min_nsnps,"

    # min mean coverage of a genomic segment
    min_cov <-",min_seg_cov,"

    # minimum value of beta above which the two alleles are present in the same number (cn equal to 1+1 or 2+2 or 3+3 ...) to account ref map bias
    equalCN.betaThr <-",equal_betaThr,"

    # Homozygous deletions threshold (change only if you know what you are doing)
    maxHomoDels <-",max_homo_dels,"

    ## Parameters of a valid deletion used to compute Adm.global
    # log2 thresholds
    deletionsLog2Levels <-",deparse(del_log_thr),"

    ## Percentage of used deletions to compute Adm.global varibility interval (change only if you know what you are doing)
    alphaPar <-",alpha_par,"

    ## Threshold on clonality value (change only if you know what you are doing)
    clonalityThreshold <-",clonal_thr,"

    ## Threshold on beta value (change only if you know what you are doing)
    betaThreshold <-",beta_thr,"

    ## Stages to perform
    # 1. analyse single sample and produe RData
    # 2. create beta table aggragating single samples analysis
    # 3. compute ploidy shift
    # 4. compute global admixture
    # 5. compute clonality table
    # 6. compute allele specific copy number table
    stages <-",deparse(stages),"

    ## CLONET can try to compute reference mapping bias and to adjust beta estimation
    computeRefMapBias <-",comp_ref_map_bias,"

    ## number of significant digits for the beta table (NA for reporting all digits)
    n.digits <-",beta_decimals,"

    ## maximum distance from (cnA, cnB) integer copy number value
    ## Dafault 0.5 corresponds to round cnA and cnB
    AllelicImbalanceTh <-",ale_imb_thr)
    cat(config_data,file=paste0(out_file_dir,"/",patient_id,".Config.R"))
  }

#' This function generates a sample_info file for CLONET pipeline
#'
#' This function takes the path to the directory with processed SNPS generated with format_SNP_data function
#' and generates a tab separated file with samples info ready to use for CLONET pipeline.
#'
#' @param patient_id Patient ID. Default Patient
#' @param snp_dir Path to formated SNP data dir.
#' @param output_dir Path to output directory.
#' @export

generate_CLONET_sample_info=function(snp_dir="",patient_id="",output_dir=""){
    sep="/"
    if(output_dir==""){
      sep=""
    }
    out_file_dir=paste0(output_dir,sep,patient_id,"_CLONET_SAMPLE_INFO")

    if (!dir.exists(out_file_dir)){
        dir.create(out_file_dir,recursive=TRUE)
    }

    out_file=paste0(patient_id,".sample.info.txt")
    files=list.files(snp_dir,full.names=TRUE,recursive=TRUE)
    tumor=files[grepl("TumorPileup",files)]
    germ=files[grepl("NormalPileup",files)]
    tumor=unlist(lapply(tumor,FUN=ULPwgs::get_sample_name))
    germ=unlist(lapply(germ,FUN=ULPwgs::get_sample_name))
    sample_info=data.frame(Tumor.Array.Name="",Tumor.Bam.Name=tumor,Normal.Array.Name="",Normal.Bam.Name=germ)
    write.table(sample_info,file=paste0(out_file_dir,"/",out_file),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
}

#' This function generates a plot of ploidy and celularity levels from CLONET data
#'
#' This function takes the path to the directory with CLONET output, as well as a tab separated file with
#' sample info and generates a plot of celularity and ploidy levels in the samples
#'
#' @param clonet_dir Path to clonet output directory
#' @param sample_data Path to file with sample info
#' @param output_dir Path to output directory.
#' @export
#' @import patchwork
#' @import tidyverse
#' @import ggplot2

plot_celullarity=function(clonet_dir="",sample_data="",output_dir=""){

    sep="/"
    if(output_dir==""){
      sep=""
    }

    if (!dir.exists(output_dir)){
        dir.create(output_dir,recursive=TRUE)
    }

    admixture=read.table(paste0(clonet_dir,"/globalAdmTable.txt"),header=TRUE,stringsAsFactors=FALSE)
    ploidy=read.table(paste0(clonet_dir,"/ploidyTable.txt"),header=TRUE,stringsAsFactors=FALSE)
    sample_info=read.table(sample_data,header=TRUE,stringsAsFactors=FALSE)
    full_data=dplyr::left_join(admixture,ploidy,by="sample")
    full_data=fuzzyjoin::fuzzy_inner_join(full_data,sample_info, by = c("sample" = "Sample_name_corrected"), match_fun = stringr::str_detect)

    plasma=full_data %>% dplyr::filter(Origin=="Plasma") %>% dplyr::mutate(Timepoint_ID=as.Date(lubridate::dmy(Timepoint_ID)))

    p=ggplot(plasma,aes(x=Timepoint_ID,y=ploidy))+geom_hline(aes(yintercept=2),linetype="dotted",alpha=0.5)+geom_bar(stat="identity",col="black",fill="red",alpha=0.5)+geom_point()+
    geom_line(aes(group=1),col="red")+theme_classic()+theme(axis.text.x = element_text(angle = 90),legend.position="bottom")+labs(x="Samples",y="Ploidy")

    p2=ggplot(plasma,aes(x=Timepoint_ID,y=1-adm))+geom_hline(aes(yintercept=0.5),linetype="dotted",alpha=0.5)+geom_bar(stat="identity",col="black",fill="blue",alpha=0.5)+geom_point()+
    geom_ribbon(data=plasma[!is.na(plasma$adm),],aes(ymin=1-adm.min, ymax=1-adm.max,group=1), fill="blue", alpha=0.2)+geom_line(data=plasma[!is.na(plasma$adm),],aes(group=1),col="blue")+theme_classic()+theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),legend.position="bottom")+labs(x="Samples",y="Celularity")+ylim(0,1)

    d=(p2/p)+plot_annotation(title = paste0(unique(full_data$Patient_ID)," Plasma"))
    ggsave(paste0(output_dir,sep,unique(full_data$Patient_ID),"_Celularity_Plasma.png"),d)

    tissue=full_data %>% dplyr::filter(Origin!="Plasma")

    if(dim(tissue)[1]>0){
      tissue$Anatomy=make.unique(tissue$Anatomy,sep="_")
      tissue$order=1-tissue$adm
      tissue$order=ifelse(is.na(tissue$order),0,tissue$order)
      p=ggplot(tissue,aes(x=reorder(Anatomy,order),y=ploidy))+geom_hline(aes(yintercept=2),linetype="dotted",alpha=0.5)+geom_bar(stat="identity",col="black",fill="red",alpha=0.5)+geom_point()+
      geom_line(aes(group=1),col="red")+theme_classic()+theme(axis.text.x = element_text(angle = 90),legend.position="bottom")+labs(x="Samples",y="Ploidy")

      p2=ggplot(tissue,aes(x=reorder(Anatomy,order),y=1-adm))+geom_hline(aes(yintercept=0.5),linetype="dotted",alpha=0.5)+geom_bar(stat="identity",col="black",fill="blue",alpha=0.5)+geom_point()+
      geom_ribbon(data=tissue[!is.na(tissue$adm),],aes(ymin=1-adm.min, ymax=1-adm.max,group=1), fill="blue", alpha=0.2)+geom_line(data=tissue[!is.na(tissue$adm),],aes(group=1),col="blue")+theme_classic()+theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x =
      element_blank(),legend.position="bottom")+labs(x="Samples",y="Celularity")+ylim(0,1)
      d=(p2/p)+plot_annotation(title = paste0(unique(full_data$Patient_ID)," Tissue"))
      ggsave(paste0(output_dir,sep,unique(full_data$Patient_ID),"_Celularity_Tissue.png"),d)
    }
}


#' This function generates a plot of allelic imbalance from CLONET output
#'
#' This function takes the path to the directory with CLONET output, as well as a tab separated file with
#' sample and gene info and generates a plot of celularity and ploidy levels in the samples
#'
#' @param clonet_dir Path to clonet output directory
#' @param sample_data Path to file with sample info
#' @param gene_data Path to file with gene info
#' @param output_dir Path to output directory.
#' @param jobs Number of jobs to run.
#' @param threads Number of threads.
#' @export
#' @import patchwork
#' @import tidyverse
#' @import ggplot2

plot_allelic_imbalance=function(clonet_dir="",sample_data="",output_dir="",gene_data="",jobs=1,threads=3){
    sep="/"
    if(output_dir==""){
      sep=""
    }

    if (!dir.exists(output_dir)){
        dir.create(output_dir,recursive=TRUE)
    }

    out_file_dir=paste0(output_dir,sep,"ALLELIC_IMBALANCE")

    if (!dir.exists(out_file_dir)){
        dir.create(out_file_dir,recursive=TRUE)
    }

    ale_imb_table=read.table(paste0(clonet_dir,"/allelicImbalanceTable.txt"),header=TRUE,stringsAsFactors=FALSE)
    clonality_table=read.table(paste0(clonet_dir,"/clonalityTable.txt"),header=TRUE,stringsAsFactors=FALSE)
    sample_info=read.table(sample_data,header=TRUE,stringsAsFactors=FALSE)
    gen_info=read.table(gene_data,header=TRUE)
    ale_imb_table=dplyr::left_join(ale_imb_table,clonality_table)

    ale_imb_table_complete=fuzzyjoin::fuzzy_inner_join(ale_imb_table, gen_info,
                    by=c("chr"="chr","start"="start","end"="end"),
                   match_fun=list(`==`, `<=`, `>=`)) %>% dplyr::mutate(Overlap="Complete")

    ale_imb_table_partially=fuzzyjoin::fuzzy_inner_join(ale_imb_table, gen_info,
                    by=c("chr"="chr","start"="end","end"="start"),
                    match_fun=list(`==`, `<=`, `>=`)) %>% dplyr::group_by(pcf_gene_symbol,sample)%>% dplyr::filter(dplyr::n()>1) %>%
                    dplyr::ungroup() %>% dplyr::mutate(Overlap="Partial")

    ale_imb_table=rbind(ale_imb_table_complete,ale_imb_table_partially)

    ale_imb_table=ale_imb_table %>% dplyr::mutate(Type=ifelse(cnA.int==1 & cnB.int==1,"WT",
      ifelse(cnA.int==1 & cnB.int==0,"Hem.Del",
      ifelse(cnA.int==2 & cnB.int==1,"Gain",
      ifelse(cnA.int==2 & cnB.int==2,"Amp",
      ifelse(cnA.int==2 & cnB.int==0,"CN‐LOH",
      ifelse(cnA.int==0 & cnB.int==0,"Hom.Del",
      ifelse(cnA.int==3 & cnB.int==0,"Gain‐LOH",
      ifelse(cnA.int>=4 & cnB.int==0,"Amp-LOH",
      ifelse(cnA.int>=3 & cnB.int==1,"Unbalanced Gain",
      ifelse((is.na(cnA.int) & is.na(cnB.int)),"NA","Unbalanced Amp")))))))))))
    ale_imb_table=ale_imb_table %>% dplyr::mutate(col=ifelse(cnA.int==1 & cnB.int==1,"grey",
      ifelse(cnA.int==1 & cnB.int==0,"lightblue",
      ifelse(cnA.int==2 & cnB.int==1,"red",
      ifelse(cnA.int==2 & cnB.int==2,"firebrick",
      ifelse(cnA.int==2 & cnB.int==0,"yellow",
      ifelse(cnA.int==0 & cnB.int==0,"blue",
      ifelse(cnA.int==3 & cnB.int==0,"orange",
      ifelse(cnA.int>=4 & cnB.int==0,"darkorange",
      ifelse(cnA.int>=3 & cnB.int==1,"brown4",
      ifelse((is.na(cnA.int) & is.na(cnB.int)),"black","deeppink4")))))))))))

    full_data=fuzzyjoin::fuzzy_inner_join(ale_imb_table,sample_info, by = c("sample" = "Sample_name_corrected"), match_fun = stringr::str_detect)
    full_data$Allelic_Imbalance=forcats::fct_rev(as.factor(ifelse(full_data$AllelicImbalance<=0.2,"E[AI] ≤ 0.2     ","E[AI] > 0.2     ")))
    full_data$Symbol=ifelse(grepl("CONTROL",full_data$pcf_gene_class),paste0(full_data$pcf_gene_symbol,"[C:","'",full_data$chr.y,substring(full_data$band,1,1),"'","]"),paste0(full_data$pcf_gene_symbol,"[T:","'",full_data$chr.y,substring(full_data$band,1,1),"'","]"))
    full_data$ID=ifelse(full_data$Origin=="Plasma",as.character(lubridate::dmy(full_data$Timepoint_ID)),full_data$Anatomy)
    full_data$AI=ifelse(full_data$AllelicImbalance<=0.2,"-",paste0(full_data$cnA.int,"/",full_data$cnA.int))
    write.table(file=paste0(out_file_dir,"/",unique(full_data$Patient_ID),".CLONET.txt"),x=full_data,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

    min.log2=min(full_data$log2,na.rm=TRUE)
    max.log2=max(full_data$log2,na.rm=TRUE)
    min.beta=0
    max.beta=1
    max.cnA=max(full_data$cnA.int,na.rm=TRUE)

  parallel::mclapply(unique(full_data$Origin),FUN=function(y){
    sub_origin=full_data %>% dplyr::filter(Origin==y)
    parallel::mclapply(unique(sub_origin$ID),FUN=function(x){
      sub_ID=sub_origin %>% dplyr::filter(ID==x)
      p01=ggplot(sub_ID,aes(x=log2))+geom_histogram(aes(y=..density..),binwidth=0.1,alpha=0.9,col="black")+geom_density(aes(y=..density..))+scale_fill_identity()+theme_classic()+xlim(min.log2-0.1,max.log2+0.1)
      pa1=ggplot(sub_ID)+geom_point(aes(y=beta,x=log2,col=col))+theme_classic()+geom_vline(aes(xintercept=0),linetype="dashed")+geom_vline(aes(xintercept=1),linetype="dashed")+geom_vline(aes(xintercept=-1),linetype="dashed")+geom_vline(aes(xintercept=0.6),linetype="dashed")+scale_color_identity()+xlim(min.log2-0.1,max.log2+0.1)+ylim(min.beta-0.1,max.beta+0.1)
      p1a=p01/pa1+plot_layout(height=c(2,8))
      p3=ggplot(sub_ID)+geom_vline(aes(xintercept=0),linetype="dashed")+geom_vline(aes(xintercept=1),linetype="dashed")+geom_vline(aes(xintercept=-1),linetype="dashed")+geom_vline(aes(xintercept=0.6),linetype="dashed")+ggrepel::geom_label_repel(data=sub_ID %>% dplyr::filter(AllelicImbalance>0.2,!is.na(pcf_gene_class)),aes(y=beta,x=log2,col=col,label=Symbol),force=20,max.overlaps=1000,min.segment.length = 0,parse=TRUE)+geom_point(data=sub_ID %>% dplyr::filter(!is.na(pcf_gene_symbol),!is.na(Allelic_Imbalance)),aes(shape=Allelic_Imbalance,y=beta,x=log2,col=col))+theme_classic()+scale_color_identity()+xlim(min.log2-0.1,max.log2+0.1)+ylim(min.beta-0.1,max.beta+0.1)+theme(legend.position=c(0.88,0.075),legend.background = element_rect(
                                        size=0.2, linetype="solid",colour="black"),
              legend.key.width=unit(1.1,"cm"))
      dummy=ggplot(sub_origin[!is.na(sub_origin$cnA),])+geom_hline(aes(yintercept=0),linetype="dashed")+geom_hline(aes(yintercept=1),linetype="dashed")+geom_hline(aes(yintercept=2),linetype="dashed")+geom_hline(aes(yintercept=3),linetype="dashed") +geom_vline(aes(xintercept=0),linetype="dashed")+geom_vline(aes(xintercept=1),linetype="dashed")+geom_vline(aes(xintercept=2),linetype="dashed")+geom_vline(aes(xintercept=3),linetype="dashed")+ geom_point(aes(x=cnA,y=cnB,col=col))+theme_classic()+geom_abline(intercept = 0, slope = 1)+xlim(0,3) +ylim(0,3)+ scale_color_identity(name = "Copy Number",
                                breaks = sub_origin$col,
                                labels = sub_origin$Type,
                                guide = "legend")
      leg=cowplot::get_legend(dummy)
      p1=ggplot(sub_ID %>% dplyr::filter(!is.na(cnA)),aes(x=cnA,y=cnB,col=col))+ geom_point()+scale_color_identity()+geom_abline(intercept = 0, slope = 1)+geom_hline(yintercept=c(1:(max.cnA)),linetype="dashed")+geom_vline(xintercept=c(1:(max.cnA)),linetype="dashed")+theme_classic()+xlim(0,max.cnA) +ylim(0,max.cnA)

      leg=ggplotify::as.ggplot(leg)
      p1a=p01/pa1
      p1b=leg/p1
      pa=(p1a|p1b)
      pa=(pa+p3)+plot_layout(width=c(3,3,6))+plot_annotation(title = x)
      ggsave(paste0(out_file_dir,"/",unique(sub_ID$Patient_ID),".",x,".CLONET_per_gene_",y,".png"),pa,width=20,height=10)
      write.table(file=paste0(out_file_dir,"/",unique(sub_ID$Patient_ID),".",x,".Allelic_Imbalance_CLONET_",y,".txt"),x=sub_ID,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
    },mc.cores=jobs)
  },mc.cores=threads)

    tc_and_ploidy_per_sample=full_data %>% dplyr::distinct(ID,adm,adm.max,adm.min,ploidy)
    print(tc_and_ploidy_per_sample)
    log2_corr_per_gene=full_data %>% dplyr::group_by(Symbol,ID) %>% dplyr::summarise(meanLog2corr=mean(log2.corr))
    AI_per_gene=full_data %>% dplyr::group_by(Symbol,ID) %>% dplyr::summarise(AI=paste0(AI,collapse="\n"))
    log2_corr_per_gene_wider=log2_corr_per_gene %>% tidyr::pivot_wider(id_cols="ID",names_from="Symbol",values_from="meanLog2corr")
    AI_per_gene_wider=AI_per_gene %>% tidyr::pivot_wider(id_cols="ID",names_from="Symbol",values_from="AI")
    write.table(file=paste0(out_file_dir,"/",unique(full_data$Patient_ID),".Allelic_Imbalance_CLONET_log2_corr_matrix.txt"),x=log2_corr_per_gene_wider,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
    write.table(file=paste0(out_file_dir,"/",unique(full_data$Patient_ID),".Allelic_Imbalance_CLONET_log2_corr_matrix.txt"),x=AI_per_gene_wider,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
    c_names=as.character(unique(log2_corr_per_gene$Symbol))
    log2_corr_per_gene_wider=as.data.frame(log2_corr_per_gene_wider)
    AI_per_gene_wider=as.data.frame(AI_per_gene_wider)
    r_names=log2_corr_per_gene_wider[rowSums(is.na(log2_corr_per_gene_wider))<(ncol(log2_corr_per_gene_wider)-1),1]
    log2_corr_mtx=log2_corr_per_gene_wider[rowSums(is.na(log2_corr_per_gene_wider))<(ncol(log2_corr_per_gene_wider)-1),-1]

    AI_per_gene_wider=AI_per_gene_wider[r_names,c_names]

    print(AI_per_gene_wider)
    tc_and_ploidy_per_sample=tc_and_ploidy_per_sample %>% tidyr::drop_na()
    rownames(tc_and_ploidy_per_sample)=tc_and_ploidy_per_sample$ID
    tc_and_ploidy_per_sample=tc_and_ploidy_per_sample[r_names,]
    col_fun = circlize::colorRamp2(c(0, 1), c("white", "red"))
    col_fun2 = circlize::colorRamp2(c(0,max(tc_and_ploidy_per_sample$ploidy)), c("blue", "brown"))
    rows_ha = ComplexHeatmap::rowAnnotation(ploidy = ComplexHeatmap::anno_barplot(tc_and_ploidy_per_sample$ploidy,baseline=2,col=col_fun2), tmf =ComplexHeatmap::anno_simple(1-tc_and_ploidy_per_sample$adm,col = col_fun))
    png(paste0(out_file_dir,"/",unique(full_data$Patient_ID),".log2_corrected_CLONET.png"),width=12,height=8,res=1200,units="in",type="cairo-png")
    ComplexHeatmap::draw(ComplexHeatmap::Heatmap(as.matrix(log2_corr_mtx),na_col="black",column_labels=c_names,
    row_labels=r_names,column_names_gp=grid::gpar(fontsize=7),cluster_rows=FALSE,cluster_columns=TRUE,
    row_split=c(rep("Plasma",sum(!grepl("[aA-zZ]",r_names))),rep("Tissue",sum(grepl("[aA-zZ]",r_names)))),name="log2.cor",right_annotation = rows_ha, layer_fun = function(j, i, x, y, width, height, fill) {
        grid::grid.text(sprintf("%.1f", ComplexHeatmap::pindex(AI_per_gene_wider,i, j)), x, y, gp = grid::gpar(fontsize = 4))
      }))
    dev.off()
}

#' This function generates a plot of ploidy and celularity levels from CLONET data
#'
#' This function takes the path to the directory with CLONET output, as well as a tab separated file with
#' sample info and generates a plot of celularity and ploidy levels in the samples
#'
#' @param cn_call_data Path to segment data with cn calls
#' @param sample_data Path to file with sample info
#' @param output_dir Path to output directory.
#' @export
#' @import patchwork
#' @import tidyverse
#' @import ggplot2

plot_cn_calls=function(cn_call_data="",sample_data="",output_dir=""){

    sep="/"
    if(output_dir==""){
      sep=""
    }

    if (!dir.exists(output_dir) & output_dir!=""){
        dir.create(output_dir,recursive=TRUE)
    }

    cn_info=read.table(cn_call_data,header=TRUE,stringsAsFactors=FALSE)
    sample_info=read.table(sample_data,header=TRUE,stringsAsFactors=FALSE)
    full_data=fuzzyjoin::fuzzy_inner_join(cn_info,sample_info, by = c("sample" = "Sample_name_corrected"), match_fun = stringr::str_detect)
    full_data$cn=ifelse(full_data$chr=="X",  full_data$cn+1,  full_data$cn)
    full_data=full_data %>% dplyr::mutate(CN=ifelse(cn>2,"GAIN",ifelse(cn<2,"LOSS","NEUTRAL")))
    full_data=full_data %>% dplyr::mutate(CNs=ifelse(CN=="GAIN"|CN=="LOSS","CNA","NEUTRAL"))
    full_data$ID=ifelse(full_data$Origin=="Plasma",as.character(lubridate::dmy(full_data$Timepoint_ID)),full_data$Anatomy)

    plasma=full_data %>% dplyr::filter(Origin=="Plasma")
    p1=ggplot(plasma %>% dplyr::group_by(sample,CN,ID,CNs)%>% dplyr::summarise(Count=dplyr::n()))+geom_bar(stat="identity",aes(x=CNs,y=Count,fill=CN),col="black")+facet_grid(.~ID)+
    theme_classic()+theme(axis.text.x = element_text(angle = 90))+plot_annotation(title=paste0(unique(plasma$Patient_ID)," Plasma"),subtitle="Somatic copy number uncorrected") +theme(strip.text.x = element_text(size = 6))
    ggsave(paste0(output_dir,sep,unique(plasma$Patient_ID),"_SCNA_count_Plasma.png"),p1,height=10,width=20)


    tissue=full_data %>% dplyr::filter(Origin!="Plasma")

    if(dim(tissue)[1]>0){
      p2=ggplot(tissue %>% dplyr::group_by(sample,CN,ID,CNs) %>% dplyr::summarise(Count=dplyr::n()) %>% dplyr::group_by(ID)%>% dplyr::mutate(TotalCN=sum(Count[CN!="NEUTRAL"])))+
      geom_bar(stat="identity",aes(x=CNs,y=Count,fill=CN),col="black")+facet_grid(~reorder(ID,TotalCN))+theme_classic()+theme(axis.text.x = element_text(angle = 90))+plot_annotation(title=paste0(unique(tissue$Patient_ID), " Tissue"),subtitle="Somatic copy number uncorrected")+theme(strip.text.x = element_text(size = 6))
      ggsave(paste0(output_dir,sep,unique(tissue$Patient_ID),"_SCNA_count_Tissue.png"),p2,height=10,width=20)
    }
}


#' This function generates a plot of reltionship between samples based on Hammering Distance
#' computed as the number of bins thats differ between samples

#' This function takes the path to the directory with CLONET output, as well as a tab separated file with
#' sample info and generates a plot of celularity and ploidy levels in the samples
#'
#' @param cn_call_data Path to segment data with cn calls
#' @param sample_data Path to file with sample info
#' @param ref_bins  Path to reference genome binned into 100kb/500kb bin
#' @param output_dir Path to output directory.
#' @param threads Number of threads to use.
#' @export
#' @import patchwork
#' @import tidyverse
#' @import ggplot2

plot_evolutionary_distance=function(cn_call_data="",sample_data="",ref_bins="",output_dir="",threads=3){

    sep="/"
    if(output_dir==""){
      sep=""
    }

    if (!dir.exists(output_dir)& output_dir!=""){
        dir.create(output_dir,recursive=TRUE)
    }
    cn_info=read.table(cn_call_data,header=TRUE,stringsAsFactors=FALSE)
    sample_info=read.table(sample_data,header=TRUE,stringsAsFactors=FALSE)
    segments=read.table(ref_bins,stringsAsFactors=FALSE)
    full_data=fuzzyjoin::fuzzy_inner_join(cn_info,sample_info, by = c("sample" = "Sample_name_corrected"), match_fun = stringr::str_detect)
    full_data$cn=ifelse(full_data$chr=="X",  full_data$cn+1,  full_data$cn)
    full_data=full_data %>% dplyr::mutate(CN=ifelse(cn>2,"GAIN",ifelse(cn<2,"LOSS","NEUTRAL")))
    full_data=full_data %>% dplyr::mutate(CNs=ifelse(CN=="GAIN"|CN=="LOSS","CNA","NEUTRAL"))
    full_data$ID=ifelse(full_data$Origin=="Plasma",as.character(lubridate::dmy(full_data$Timepoint_ID)),full_data$Anatomy)
    segments$cn=2
    segments$region="-"
    solution=parallel::mclapply(unique(full_data$ID),FUN=function(x){segments_tmp=segments;CN_sub=full_data %>% dplyr::filter(ID==x);
    for (y in 1:nrow(CN_sub)){segments_tmp[segments_tmp$V1==CN_sub[y,]$chr & segments_tmp$V2>=CN_sub[y,]$start & segments_tmp$V3<=CN_sub[y,]$end,"cn"]=CN_sub[y,]$cn;
    segments_tmp$sample=x;segments_tmp[segments_tmp$V1==CN_sub[y,]$chr & segments_tmp$V2<=CN_sub[y,]$end & segments_tmp$V3>=CN_sub[y,]$start,"region"]=paste0(CN_sub[y,]$chr,":",CN_sub[y,]$start,"-",CN_sub[y,]$end)};
    return(segments_tmp)},mc.cores=threads)
    solution=solution %>% dplyr::bind_rows()
    solution$change=ifelse(solution$cn==2,0,ifelse(solution$cn>2,1,-1))
    solution_wider=tidyr::pivot_wider(solution,id_cols="V4",names_from="sample",values_from="change")
    solution_matrix=solution_wider[,-1]
    write.table(file=paste0(output_dir,sep,unique(full_data$Patient_ID),".CN.Matrix.txt"),x=solution_wider,quote=FALSE,row.names=FALSE,col.names=TRUE)

    dist_matrix_all=dist(t(solution_matrix))
    NJ_data_all=phangorn::NJ(dist_matrix_all)
    NJ_tree_all=ape::ladderize(NJ_data_all)

    png(filename=paste0(output_dir,sep,unique(full_data$Patient_ID),".NJ_all.png"),units="in",width=12,height=10,res=1000)
      p=ape::plot.phylo(NJ_tree_all)
      p
    dev.off()



    dist_matrix_plasma=dist(t(solution_matrix[,!grepl("[aA-zZ]",colnames(solution_matrix))]))
    NJ_data_plasma=phangorn::NJ(dist_matrix_plasma)
    NJ_tree_plasma=ape::ladderize(NJ_data_plasma)

    png(filename=paste0(output_dir,sep,unique(full_data$Patient_ID),".NJ_all_plasma.png"),units="in",width=12,height=10,res=1000)
      p=ape::plot.phylo(NJ_tree_plasma)
      p
    dev.off()


    dist_matrix_tissue=dist(t(solution_matrix[,grepl("[aA-zZ]",colnames(solution_matrix))]))
    NJ_data_tissue=phangorn::NJ(dist_matrix_tissue)
    NJ_tree_tissue=ape::ladderize(NJ_data_tissue)

    png(filename=paste0(output_dir,sep,unique(full_data$Patient_ID),".NJ_all_tissue.png"),units="in",width=12,height=10,res=1000)
      p=ape::plot.phylo(NJ_tree_tissue)
      p
    dev.off()

    tp=solution_matrix
    rownames(tp)=solution_wider[,1]
    tp_pos=parallel::mcsapply(1:ncol(tp),FUN=function(x){
      tp[,x]!=dplyr::lag(tp[,x])
    },mc.cores=threads)
    tp_pos[1,]=FALSE
    dummy_tp=tp
    dummy_tp[!is.na(dummy_tp)]=0
    dummy_tp[ts_point]=1
    tp_all=dummy_tp[rowSums(dummy_tp)>1,]

    png(paste0(out_file_dir,"/",unique(full_data$Patient_ID),".transition_points.png"),width=12,height=8,res=1200,units="in",type="cairo-png")
    ComplexHeatmap::draw(ComplexHeatmap::Heatmap(tp_all,row_names_gp=grid::gpar(fontsize=5),clustering_distance_rows=function(m) dist(m,method="manhattan"),clustering_distance_columns=function(m) dist(m,method="manhattan")))
    dev.off()
}
