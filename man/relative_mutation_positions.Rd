% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_mutations.R
\name{relative_mutation_positions}
\alias{relative_mutation_positions}
\title{Relative mutation positions within the gene body}
\usage{
relative_mutation_positions(
  mutations,
  annotation,
  bin = FALSE,
  by = c(region = "gene_name"),
  threads = 1
)
}
\arguments{
\item{mutations}{either the tibble output from \code{find_variants}. It must have columns \code{seqnames},
\code{pos}, and a third column for specifying the gene id or gene name. The mutation must be within the gene region.}

\item{annotation}{Either path to the annotation file (GTF/GFF) or a GRanges object of the gene annotation.}

\item{bin}{logical(1): whether to bin the relative positions into 100 bins.}

\item{by}{character(1): the column name in the annotation to match with the gene annotation.
E.g. \code{c("region" = "gene_name")} to match the `region` column in the mutations with the
`gene_name` column in the annotation.}

\item{threads}{integer(1): number of threads to use.}
}
\value{
If \code{bin = FALSE}, a list of numeric vectors of relative positions of each mutation within the gene body.
If \code{bin = TRUE}, a numeric vector of length 100 of the number of mutations in each bin.
}
\description{
Given a set of mutations and gene annotation, calculate the relative position of each mutation
within the gene body they are in.
}
\examples{
outdir <- tempfile()
dir.create(outdir)
genome_fa <- system.file("extdata", "rps24.fa.gz", package = "FLAMES")
minimap2_align( # align to genome
  config = jsonlite::fromJSON(
    system.file("extdata", "config_sclr_nanopore_3end.json", package = "FLAMES")),
  fa_file = genome_fa,
  fq_in = system.file("extdata", "fastq", "demultiplexed.fq.gz", package = "FLAMES"),
  annot = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
  outdir = outdir
)
variants <- find_variants(
  bam_path = file.path(outdir, "align2genome.bam"),
  reference = genome_fa,
  annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
  min_nucleotide_depth = 4
)
positions <- 
 relative_mutation_positions(
   mutations = variants,
   annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
 )
}
