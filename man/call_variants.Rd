% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{call_variants}
\alias{call_variants}
\title{Variant calling of Somatic and Germline (MuTECT2, Platypus and Strelka) mutations and produces their subsequent annotation (VEP)}
\usage{
call_variants(
  bin_path = "tools/gatk/gatk",
  bin_path2 = "tools/bcftools/bcftools",
  bin_path3 = "tools/htslib/bgzip",
  bin_path4 = "tools/htslib/tabix",
  bin_path5 = "tools/platypus/Platypus.py",
  bin_path6 = "tools/ensembl-vep/vep",
  bin_path7 = "tools/strelka-2.9.10/build/bin/configureStrelkaGermlineWorkflow.py",
  bin_path8 = "tools/manta-1.6.0/build/bin/configManta.py",
  bin_path9 = "tools/svaba/bin/svaba",
  bam_dir = "",
  patient_id = "",
  germ_pattern = "GL",
  ref_genome = "",
  germ_resource = "",
  pon = "",
  output_dir = "",
  region_bed = "",
  chr_filter = "canonical",
  db = "",
  interval = "",
  targets = "",
  orientation = TRUE,
  resources = "",
  info_key = "CNN_2D",
  snp_tranche = 99.95,
  indel_tranche = 99.4,
  threads = 3,
  verbose = FALSE,
  targeted = TRUE,
  validate_id = TRUE
)
}
\arguments{
\item{bin_path}{\link{REQUIRED} Path to gatk binary. Default tools/gatk/gatk.}

\item{bin_path2}{\link{REQUIRED} Path to bcftools binary. Default tools/bcftools/bcftools.}

\item{bin_path3}{\link{REQUIRED} Path to bgzip binary. Default tools/htslib/bgzip.}

\item{bin_path4}{\link{REQUIRED} Path to tabix binary. Default tools/htslib/tabix.}

\item{bin_path5}{\link{REQUIRED} Path to platypus binary. Default tools/platypus/Platypus.py.}

\item{bin_path6}{\link{REQUIRED} Path to vep binary binary. tools/ensembl-vep/vep}

\item{bin_path7}{\link{REQUIRED} Path to Strelka Germline Workflow binary. tools/strelka-2.9.10/build/bin/configureStrelkaGermlineWorkflow.py}

\item{bin_path8}{\link{REQUIRED} Path to Manta binary. tools/manta-1.6.0/build/bin/configManta.py}

\item{bin_path9}{\link{REQUIRED} Path to svaba binary. tools/svaba/bin/svaba}

\item{bam_dir}{\link{REQUIRED} Path to directory with BAM files.}

\item{patient_id}{\link{REQUIRED} Patient ID to analyze. Has to be in file names to subselect samples.}

\item{germ_pattern}{\link{REQUIRED} Pattern used to identify germline samples. Ex GL}

\item{ref_genome}{\link{REQUIRED} Path to reference genome fasta file.}

\item{germ_resource}{\link{REQUIRED} Path to germline resources vcf file.}

\item{pon}{\link{OPTIONAL} Path to panel of normal.}

\item{output_dir}{\link{OPTIONAL} Path to the output directory.}

\item{region_bed}{\link{REQUIRED} Path to bed file with regions to analyze.}

\item{chr_filter}{\link{OPTIONAL} Chromosomes to analyze. canonical/autosomal/all or a list of chromosomes}

\item{db}{\link{OPTIONAL} Path to vcf with common variants. Used for contamination estimation.}

\item{interval}{\link{OPTIONAL} Path to interval for common variants to analyze. Used for contamination estimation.}

\item{targets}{\link{OPTIONAL} Path to BED file with capture targets.}

\item{orientation}{\link{OPTIONAL} Generate a read orientation model to filter variants. Default TRUE}

\item{info_key}{\link{OPTIONAL} Annotation column to select. Default CNN_D2}

\item{snp_tranche}{\link{OPTIONAL} SNP tranche filter value. Default 99.95}

\item{indel_tranche}{\link{OPTIONAL} Indel tranche filter value. Default 99.44}

\item{threads}{\link{OPTIONAL} Number of threads. Default 3}

\item{verbose}{\link{OPTIONAL} Enables progress messages. Default False.}

\item{targeted}{\link{OPTIONAL} Sequencing method Exome/Targeted. Default TRUE}

\item{validate_id}{\link{OPTIONAL} Validate patient_id in BAM dir. Default TRUE.}

\item{output_name}{\link{OPTIONAL} Name for the output. If not given the name of one of the samples will be used.}
}
\description{
This function calls somatic and germline mutations, in previously pre-processed sequencing data.
This function takes a path to a directory with BAM files and a patients ID. Then it subsets all BAM
files specific to the patient, identifying between cancer and normal Samples using the germline identifier.
Afterwards, it calls the MuTECT2 somatic variant caller, to indentify somatics variants, and Platypus genotyper all variants found in the samples.
Finally, it calls Varant Effect Predictor (VEP) and annotates all the variants.
}
