% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{FilterVariantTranches}
\alias{FilterVariantTranches}
\title{Filter Variant Tranches (Mutect2)}
\usage{
FilterVariantTranches(
  bin_path = "tools/gatk/gatk",
  vcf = "",
  resources = "",
  output_name = "",
  info_key = "CNN_1D",
  snp_tranche = 99.95,
  indel_tranche = 99.4,
  output_dir = "",
  keep_previous_filters = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{bin_path}{\link{REQUIRED} Path to GATK binary. Default tools/gatk/gatk}

\item{vcf}{\link{REQUIRED} Path to annotated VCF with CNN_D1/CNN_D2 scores}

\item{resources}{\link{OPTIONAL} Path to resources for variant filtering}

\item{output_name}{\link{OPTIONAL} Name of the sample to output}

\item{info_key}{\link{OPTIONAL} Annotation column to select. Default CNN_D1}

\item{snp_tranche}{\link{OPTIONAL} SNP tranche filter value. Default 99.95}

\item{indel_tranche}{\link{OPTIONAL} Indel tranche filter value. Default 99.4}

\item{output_dir}{\link{OPTIONAL} Path to output dir}

\item{keep_previous_filters}{\link{OPTIONAL} Keep previous filters in VCF. Default False}

\item{verbose}{\link{Optional} Enables progress messages. Default False}
}
\description{
This function takes an annotated VCF with CNN_D1 or CNN_D2 scores and filter the
variants based on set threshold.
}
