% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{call_sv_manta_parallel}
\alias{call_sv_manta_parallel}
\title{Call Somatic/Germline variants using MANTA in parallel
This function takes a pair of matched samples or single | multiple germline variants
and calls variants on them.
This function parallelizes the number of jobs that can be run at single time.}
\usage{
call_sv_manta_parallel(
  bin_path = "tools/manta-1.6.0/build/bin/configManta.py",
  bam_dir = "",
  ref_genome = "",
  output_dir = "",
  verbose = FALSE,
  targeted = FALSE,
  jobs = 1,
  threads = 3
)
}
\arguments{
\item{bin_path}{\link{REQUIRED} Path to strelka binary. Somatic or Germline.}

\item{bam_dir}{\link{REQUIRED} Path to directory with BAM files.}

\item{ref_genome}{\link{REQUIRED} Path to reference genome.}

\item{output_dir}{\link{OPTIONAL} Path to the output directory.}

\item{verbose}{\link{DEFAULT==FALSE} Enables progress messages.}

\item{targeted}{\link{OPTIONAL} If exome/capture method. Default FALSE}

\item{jobs}{\link{OPTIONAL} Number of jobs. Default 1}

\item{threads}{\link{OPTIONAL} Number of threads per job. Default 3}
}
\description{
Call Somatic/Germline variants using MANTA in parallel
This function takes a pair of matched samples or single | multiple germline variants
and calls variants on them.
This function parallelizes the number of jobs that can be run at single time.
}
