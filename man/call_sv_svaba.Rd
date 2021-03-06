% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{call_sv_svaba}
\alias{call_sv_svaba}
\title{Structural variant calling using svaba}
\usage{
call_sv_svaba(
  tumor_bam = "",
  bin_path = "tools/svaba/bin/svaba",
  bin_path2 = "tools/htslib/bgzip",
  bin_path3 = "tools/htslib/tabix",
  normal_bam = "",
  ref_genome = "",
  threads = 3,
  output_name = "",
  targets = "",
  dbsnp_indels = "",
  verbose = FALSE,
  output_dir = ""
)
}
\arguments{
\item{tumor_bam}{\link{REQUIRED}  Path to tumor bam file.}

\item{bin_path}{\link{REQUIRED} Path to svaba binary executable. Default path tools/svaba/svaba.}

\item{bin_path2}{\link{REQUIRED} Path to bgzip binary executable. Default path tools/htslib/bgzip.}

\item{bin_path3}{\link{REQUIRED} Path to tabix binary executable. Default path tools/htslib/tabix.}

\item{normal_bam}{\link{OPTIONAL} Path to germline bam file.}

\item{ref_genome}{\link{REQUIRED} Path to reference genome fasta file.}

\item{threads}{\link{OPTIONAL} Number of threads to use. Default 3.}

\item{output_name}{\link{OPTIONAL} Name for the output.}

\item{targets}{\link{OPTIONAL} BED file with capture target regions.}

\item{dbsnp_indels}{\link{OPTIONAL} Database with indel annotations.}

\item{verbose}{\link{OPTIONAL} Enables progress messages. Default False.}

\item{output_dir}{\link{OPTIONAL} Path to the output directory.}
}
\description{
This function calls structural variants in a pair of tumor-normal or in single tumor/normal samples matched samples
using svaba
}
