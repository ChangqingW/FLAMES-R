% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-MultiSampleSCPipeline.R
\name{MultiSampleSCPipeline}
\alias{MultiSampleSCPipeline}
\title{Pipeline for multi-sample long-read scRNA-seq data}
\usage{
MultiSampleSCPipeline(
  config_file,
  outdir,
  fastq,
  annotation,
  genome_fa,
  genome_mmi,
  minimap2,
  samtools,
  barcodes_file,
  expect_cell_number,
  controllers
)
}
\arguments{
\item{config_file}{Path to the JSON configuration file. See \code{\link{create_config}} for creating one.}

\item{outdir}{Path to the output directory. If it does not exist, it will be created.}

\item{fastq}{A named vector of fastq file (or folder) paths. Each element of the vector
will be treated as a sample. The names of the vector will be used as the sample names.
If not named, the sample names will be generated from the file names.}

\item{annotation}{The file path to the annotation file in GFF3 / GTF format.}

\item{genome_fa}{The file path to the reference genome in FASTA format.}

\item{genome_mmi}{(optional) The file path to minimap2's index reference genome.}

\item{minimap2}{(optional) The path to the minimap2 binary. If not provided, FLAMES will
use a copy from bioconda via \code{basilisk}.}

\item{samtools}{(optional) The path to the samtools binary. If not provided, FLAMES will
use a copy from bioconda via \code{basilisk}.}

\item{barcodes_file}{The file with expected cell barcodes, with each barcode on a new line.}

\item{expect_cell_number}{The expected number of cells in the sample. This is used if
\code{barcodes_file} is not provided. See \code{BLAZE} for more details.}

\item{controllers}{(optional, \strong{experimental}) A \code{crew_class_controller} object for running certain steps}
}
\value{
A \code{FLAMES.MultiSampleSCPipeline} object. The pipeline can be run using
the \code{\link{run_FLAMES}} function. The resulting list of SingleCellExperiment
objects can be accessed using the \code{experiment} method.
}
\description{
Semi-supervised isofrom detection and annotation for long read data. This variant is
meant for multi-sample scRNA-seq data. Specific parameters can be configured in
the config file (see \code{\link{create_config}}), input files are specified via
arguments.
}
\details{
By default the pipeline starts with demultiplexing the input fastq data. If the
cell barcodes are known apriori (e.g. via coupled short-read sequencing), the
\code{barcodes_file} argument can be used to specify a file containing the cell
barcodes, and a modified Rcpp version of \code{flexiplex} will be used; otherwise,
\code{expect_cell_number} need to be provided, and \code{BLAZE} will be used to
generate the cell barcodes. The pipeline then aligns the reads to the genome using
\code{minimap2}. The alignment is then used for isoform detection (either using
\code{FLAMES} or \code{bambu}, can be configured). The reads are then realigned
to the detected isoforms. Finally, a transcript count matrix is generated (either
using \code{FLAMES}'s simplistic counting or \code{oarfish}'s Expectation
Maximization algorithm, can be configured). The results can be accssed with
\code{experiment(pipeline)}. If the pipeline errored out / new steps were configured,
it can be resumed by calling \code{resume_FLAMES(pipeline)}
}
\examples{
reads <- ShortRead::readFastq(
  system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES")
)
outdir <- tempfile()
dir.create(outdir)
dir.create(file.path(outdir, "fastq"))
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
ShortRead::writeFastq(reads[1:100],
  file.path(outdir, "fastq/sample1.fq.gz"), mode = "w", full = FALSE)
reads <- reads[-(1:100)]
ShortRead::writeFastq(reads[1:100],
  file.path(outdir, "fastq/sample2.fq.gz"), mode = "w", full = FALSE)
reads <- reads[-(1:100)]
ShortRead::writeFastq(reads,
  file.path(outdir, "fastq/sample3.fq.gz"), mode = "w", full = FALSE)
ppl <- MultiSampleSCPipeline(
  config_file = create_config(outdir, type = "sc_3end", threads = 1, no_flank = TRUE),
  outdir = outdir,
  fastq = c("sampleA" = file.path(outdir, "fastq"),
    "sample1" = file.path(outdir, "fastq", "sample1.fq.gz"),
    "sample2" = file.path(outdir, "fastq", "sample2.fq.gz"),
    "sample3" = file.path(outdir, "fastq", "sample3.fq.gz")),
  annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
  genome_fa = genome_fa,
  barcodes_file = rep(bc_allow, 4)
)
ppl <- run_FLAMES(ppl)
experiment(ppl)
}
\seealso{
\code{\link{SingleCellPipeline}} for single-sample long data and more details on the
pipeline output,
\code{\link{create_config}} for creating a configuration file,
\code{\link{BulkPipeline}} for bulk long data.
}
