% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{filter_VEP}
\alias{filter_VEP}
\title{Filter VEP results}
\usage{
filter_VEP(
  bin_path = "tools/ensembl-vep/filter_vep",
  bin_path2 = "tools/htslib/bgzip",
  bin_path3 = "tools/htslib/tabix",
  unf_vcf = "",
  filter = "",
  verbose = FALSE,
  output_dir = ""
)
}
\arguments{
\item{bin_path}{Path to vep binary. Default tools/ensembl-vep/filter_vep.}

\item{bin_path2}{Path to bgzip binary. Default tools/htslib/bgzip.}

\item{bin_path3}{Path to tabix binary. Default tools/htslib/tabix.}

\item{unf_vcf}{Path to vcf file.}

\item{filter}{Column tags to filter by. Examples "(MAX_AF < 0.01 or not MAX_AF) and FILTER = PASS and SYMBOL in /home/regmova/PCFselect/Panel_V2/PANEL_GENE.txt". For more information check VEP filter info.}

\item{verbose}{Enables progress messages. Default False.}

\item{output_dir}{Path to the output directory.}
}
\description{
This function filters VCFs generated by VEP
}
