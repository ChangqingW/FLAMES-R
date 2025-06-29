% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_decay.R
\name{get_coverage}
\alias{get_coverage}
\title{Get read coverages from BAM file}
\usage{
get_coverage(bam, min_counts = 10, remove_UTR = FALSE, annotation)
}
\arguments{
\item{bam}{path to the BAM file, or a parsed GAlignments object}

\item{min_counts}{numeric, the minimum number of alignments required for a transcript to be included}

\item{remove_UTR}{logical, if \code{TRUE}, remove the UTRs from the coverage}

\item{annotation}{(Required if \code{remove_UTR = TRUE}) path to the GTF annotation file}
}
\value{
a tibble of the transcript information and coverages, with the following columns:
\itemize{
\item \code{transcript}: the transcript name / ID
\item \code{read_counts}: the total number of aligments for the transcript
\item \code{coverage_1-100}: the coverage at each of the 100 positions along the transcript
\item \code{tr_length}: the length of the transcript
}
}
\description{
Get the read coverages for each transcript in the BAM file (or a GAlignments object).
The read coverages are sampled at 100 positions along the transcript, and the coverage is scaled
by dividing the coverage at each position by the total read counts for the transcript. If a BAM
file is provided, alignment with MAPQ < 5, secondary alignments and supplementary alignments
are filtered out. A \code{GAlignments} object can also be provided in case alternative filtering
is desired.
}
\examples{
ppl <- example_pipeline("BulkPipeline")
steps(ppl)["isoform_identification"] <- FALSE
ppl <- run_step(ppl, "read_realignment")
x <- get_coverage(ppl@transcriptome_bam[[1]])
head(x)
}
