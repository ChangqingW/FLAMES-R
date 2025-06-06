% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BLAZE_demultiplexing.R
\name{blaze}
\alias{blaze}
\title{BLAZE Assign reads to cell barcodes.}
\usage{
blaze(expect_cells, fq_in, ...)
}
\arguments{
\item{expect_cells}{Integer, expected number of cells. Note: this could be just a rough estimate. E.g., the targeted number of cells.}

\item{fq_in}{File path to the fastq file used as a query sequence file}

\item{...}{Additional BLAZE configuration parameters. E.g., setting
`'output-prefix'='some_prefix'` is equivalent to specifying `--output-prefix some_prefix` in BLAZE; Similarly,
`overwrite=TRUE` is equivalent to switch on the `--overwrite` option. Note that the specified parameters will
override the parameters specified in the configuration file. All available options can be found at https://github.com/shimlab/BLAZE.}
}
\value{
A \code{data.frame} summarising the reads aligned. Other outputs are written to disk. 
The details of the output files can be found at https://github.com/shimlab/BLAZE.
}
\description{
Uses BLAZE to generate barcode list and assign reads to cell barcodes.
}
\examples{
outdir <- tempfile()
dir.create(outdir)
fastq <- system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES")
blaze(
  expect_cells = 10, fastq,
  "output-prefix" = file.path(outdir, ""),
  "output-fastq" = file.path(outdir, "output.fastq"),
  overwrite=TRUE
)

}
