% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{vcf_filtering}
\alias{vcf_filtering}
\title{VCF filtering using GATK}
\usage{
vcf_filtering(
  bin_path = "tools/gatk/gatk",
  bin_path2 = "tools/htslib/bgzip",
  bin_path3 = "tools/htslib/tabix",
  unfil_vcf = "",
  ref_genome = "",
  unfil_vcf_stats = "",
  contamination = "",
  segmentation = "",
  orientation_model = "",
  verbose = FALSE,
  output_dir = ""
)
}
\arguments{
\item{bin_path}{\link{REQUIRED} Path to gatk binary. Default tools/gatk/gatk.}

\item{bin_path2}{\link{REQUIRED} Path to bgzip binary. Default tools/htslib/bgzip.}

\item{bin_path3}{\link{REQUIRED} Path to tabix binary. Default tools/htslib/tabix.}

\item{unfil_vcf}{\link{REQUIRED} Path to unfiltered vcf file.}

\item{ref_genome}{\link{REQUIRED} Path to reference genome fasta file.}

\item{unfil_vcf_stats}{\link{REQUIRED} Path to unfiltered vcf file stats.}

\item{contamination}{\link{OPTIONAL} Path to contamination table. Also requires segmentation.}

\item{segmentation}{\link{OPTIONAL} Path to segments table. Also requires contamination.}

\item{orientation_model}{\link{OPTIONAL} Path to orientation model generated using learn_orientation.}

\item{verbose}{\link{OPTIONAL} Enables progress messages. Default False.}

\item{output_dir}{\link{OPTIONAL} Path to the output directory.}
}
\description{
This function filters VCF calls using GATK statistics
}
