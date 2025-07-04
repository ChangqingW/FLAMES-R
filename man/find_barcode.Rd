% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_barcode.R
\name{find_barcode}
\alias{find_barcode}
\title{Match Cell Barcodes}
\usage{
find_barcode(
  fastq,
  barcodes_file,
  max_bc_editdistance = 2,
  max_flank_editdistance = 8,
  reads_out,
  stats_out,
  threads = 1,
  pattern = c(primer = "CTACACGACGCTCTTCCGATCT", BC = paste0(rep("N", 16), collapse =
    ""), UMI = paste0(rep("N", 12), collapse = ""), polyT = paste0(rep("T", 9), collapse
    = "")),
  TSO_seq = "",
  TSO_prime = 3,
  strand = "+",
  cutadapt_minimum_length = 1,
  full_length_only = FALSE
)
}
\arguments{
\item{fastq}{A path to a FASTQ file or a directory containing FASTQ files.}

\item{barcodes_file}{path to file containing barcode allow-list, with one barcode in each line}

\item{max_bc_editdistance}{max edit distances for the barcode sequence}

\item{max_flank_editdistance}{max edit distances for the flanking sequences (primer and polyT)}

\item{reads_out}{path to output FASTQ file}

\item{stats_out}{path of output stats file}

\item{threads}{number of threads to be used}

\item{pattern}{named character vector defining the barcode pattern}

\item{TSO_seq}{TSO sequence to be trimmed}

\item{TSO_prime}{either 3 (when \code{TSO_seq} is on 3' the end) or 5 (on 5' end)}

\item{strand}{strand of the barcode pattern, either '+' or '-' (read will be reverse complemented
after barcode matching if '-')}

\item{cutadapt_minimum_length}{minimum read length after TSO trimming (cutadapt's --minimum-length)}

\item{full_length_only}{boolean, when TSO sequence is provided, whether reads without TSO
are to be discarded}
}
\value{
a list containing: \code{reads_tb} (tibble of read demultiplexed information) and
\code{input}, \code{output}, \code{read1_with_adapter} from cutadapt report
(if TSO trimming is performed)
}
\description{
demultiplex reads with flexiplex
}
\details{
This function demultiplexes reads by searching for flanking sequences (adaptors)
around the barcode sequence, and then matching against allowed barcodes.
}
\examples{
outdir <- tempfile()
dir.create(outdir)
bc_allow <- file.path(outdir, "bc_allow.tsv")
R.utils::gunzip(
  filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
  destname = bc_allow, remove = FALSE
)
find_barcode(
  fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
  stats_out = file.path(outdir, "bc_stat.tsv.gz"),
  reads_out = file.path(outdir, "demultiplexed.fastq.gz"),
  barcodes_file = bc_allow, 
  TSO_seq = "AAGCAGTGGTATCAACGCAGAGTACATGGG", TSO_prime = 5,
  strand = '-', cutadapt_minimum_length = 10, full_length_only = TRUE
)
}
