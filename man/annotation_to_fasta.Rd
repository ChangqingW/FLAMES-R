% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_isoform.R
\name{annotation_to_fasta}
\alias{annotation_to_fasta}
\title{GTF/GFF to FASTA conversion}
\usage{
annotation_to_fasta(isoform_annotation, genome_fa, outfile, extract_fn)
}
\arguments{
\item{isoform_annotation}{Path to the annotation file (GTF/GFF3)}

\item{genome_fa}{The file path to genome fasta file.}

\item{outfile}{The file path to the output FASTA file.}

\item{extract_fn}{(optional) Function to extract a \code{GRangesList} object
E.g. \code{function(grl){GenomicFeatures::cdsBy(grl, by="tx")}}}
}
\value{
This does not return anything. A FASTA file will be created at the specified location.
}
\description{
convert the transcript annotation to transcriptome assembly as FASTA file.
}
\examples{
fasta <- tempfile()
annotation_to_fasta(system.file("extdata", "rps24.gtf.gz", package = "FLAMES"), system.file("extdata", "rps24.fa.gz", package = "FLAMES"), fasta)
cat(readChar(fasta, 1e3))

}
