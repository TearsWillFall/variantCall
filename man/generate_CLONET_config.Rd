% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{generate_CLONET_config}
\alias{generate_CLONET_config}
\title{This function generates a config files for CLONET pipeline}
\usage{
generate_CLONET_config(
  patient_id = "PATIENT",
  clonet_dir = "tools/CLONET",
  snp_dir = "",
  segment_data = "",
  sample_info = "",
  min_snp_cov = 10,
  min_nsnps = 10,
  min_seg_cov = 20,
  equal_betaThr = 0.9,
  max_homo_dels = 0.01,
  del_log_thr = c(-1, -0.25),
  alpha_par = 0.9,
  clonal_thr = 0.85,
  beta_thr = 0.85,
  stages = c(1, 2, 3, 4, 5, 6),
  comp_ref_map_bias = FALSE,
  beta_decimals = 3,
  ale_imb_thr = 0.5,
  beta_method = "STM",
  adm_method = "2D",
  jobs = 1,
  threads = 3,
  output_dir = ""
)
}
\arguments{
\item{patient_id}{Patient ID. Default Patient}

\item{clonet_dir}{Path to CLONET dir}

\item{snp_dir}{Path to formated SNP data dir.}

\item{segment_data}{Path to file with formated segment data.}

\item{sample_info}{Path to file with sample info data.}

\item{min_snp_cov}{Minimum tumor coverage for informative SNPs. Default 10.}

\item{min_nsnps}{Minimum number of SNPs per segment. Default 10.}

\item{min_seg_cov}{Minimum segment coverage. Default 20.}

\item{equal_betaThr}{Minimum value of beta above which the two alleles are present in the same number. Default 0.9}

\item{max_homo_dels}{Homozygous deletions threshold. Default 0.01}

\item{del_log_thr}{Parameters of a valid deletion used to compute Adm.global. Default c(-1,-0.25)}

\item{alpha_par}{Percentage of used deletions to compute Adm.global varibility interval. Default 0.9}

\item{clonal_thr}{Clonality value threshold. Default 0.85}

\item{beta_thr}{Beta value threshold. Default 0.85}

\item{stages}{Analysis stages to run trough.Default c(1,2,3,4,5,6)}

\item{comp_ref_map_bias}{Compute reference mapping bias and to adjust beta estimation. Default FALSE}

\item{beta_decimals}{Number of beta value decimals to report. Default 3}

\item{beta_method}{Method for beta estimation. Default STM. Options STM/GB}

\item{adm_method}{Method for admixture estimation. Default 2D. Options 1D/2D}

\item{jobs}{Number of Samples to analyze in parallel. Default 1}

\item{threads}{Number of threads per job to use. Default 3}

\item{output_dir}{Path to the output directory.}

\item{verbose}{Enables progress messages. Default FALSE.}
}
\description{
This function takes multiple parameters and generates a config file for the CLONET tool
}
