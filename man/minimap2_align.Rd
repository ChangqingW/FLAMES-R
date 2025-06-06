% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minimap2_align.R
\name{minimap2_align}
\alias{minimap2_align}
\title{Minimap2 Align to Genome}
\usage{
minimap2_align(
  fq_in,
  fa_file,
  config,
  outfile,
  minimap2_args,
  sort_by,
  minimap2,
  samtools,
  threads = 1,
  tmpdir
)
}
\arguments{
\item{fq_in}{File path to the fastq file used as a query sequence file}

\item{fa_file}{Path to the fasta file used as a reference database for alignment}

\item{config}{Parsed list of FLAMES config file}

\item{outfile}{Path to the output file}

\item{minimap2_args}{Arguments to pass to minimap2, see minimap2 documentation for details.}

\item{sort_by}{Column to sort the bam file by, see \code{samtools sort} for details}

\item{minimap2}{Path to minimap2 binary}

\item{samtools}{path to the samtools binary.}

\item{threads}{Integer, threads for minimap2 to use, see minimap2 documentation for details,}

\item{tmpdir}{Temporary directory to use for intermediate files.
FLAMES will try to detect cores if this parameter is not provided.}
}
\value{
a \code{data.frame} summarising the reads aligned
}
\description{
Uses minimap2 to align sequences agains a reference databse.
Uses options '-ax splice -t 12 -k14 --secondary=no \code{fa_file} \code{fq_in}'
}
\keyword{internal}
