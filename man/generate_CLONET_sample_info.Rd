% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{generate_CLONET_sample_info}
\alias{generate_CLONET_sample_info}
\title{This function generates a sample_info file for CLONET pipeline}
\usage{
generate_CLONET_sample_info(snp_dir = "", patient_id = "", output_dir = "")
}
\arguments{
\item{snp_dir}{Path to formated SNP data dir.}

\item{patient_id}{Patient ID. Default Patient}

\item{output_dir}{Path to output directory.}
}
\description{
This function takes the path to the directory with processed SNPS generated with format_SNP_data function
and generates a tab separated file with samples info ready to use for CLONET pipeline.
}
