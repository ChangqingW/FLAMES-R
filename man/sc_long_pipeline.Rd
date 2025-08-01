% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-SingleCellPipeline.R
\name{sc_long_pipeline}
\alias{sc_long_pipeline}
\title{Pipeline for Single Cell Data (deprecated)}
\usage{
sc_long_pipeline(
  annotation,
  fastq,
  outdir,
  genome_fa,
  minimap2 = NULL,
  barcodes_file = NULL,
  expect_cell_number = NULL,
  config_file = NULL
)
}
\arguments{
\item{annotation}{The file path to the annotation file in GFF3 format}

\item{fastq}{The file path to input fastq file}

\item{outdir}{The path to directory to store all output files.}

\item{genome_fa}{The file path to genome fasta file.}

\item{minimap2}{Path to minimap2, optional.}

\item{barcodes_file}{The file with expected cell barcodes, with each barcode on a new line.}

\item{expect_cell_number}{The expected number of cells in the sample. This is used if
\code{barcodes_file} is not provided. See \code{BLAZE} for more details.}

\item{config_file}{File path to the JSON configuration file.}
}
\value{
A \code{SingleCellPipeline} object containing the transcript counts.
}
\description{
This function is deprecated. Please use [SingleCellPipeline()] instead.
}
\examples{
outdir <- tempfile()
dir.create(outdir)
bc_allow <- file.path(outdir, "bc_allow.tsv")
genome_fa <- file.path(outdir, "rps24.fa")
R.utils::gunzip(
  filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
  destname = bc_allow, remove = FALSE
)
R.utils::gunzip(
  filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
  destname = genome_fa, remove = FALSE
)
sce <- FLAMES::sc_long_pipeline(
  genome_fa = genome_fa,
  fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
  annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
  outdir = outdir,
  barcodes_file = bc_allow
)
}
\seealso{
\code{\link{SingleCellPipeline}} for the new pipeline interface,
\code{\link{BulkPipeline}} for bulk long data,
\code{\link{MultiSampleSCPipeline}} for multi sample single cell pipelines.
}
