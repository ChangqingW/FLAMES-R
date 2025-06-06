% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_annotate_plots.R
\name{plot_isoform_reduced_dim}
\alias{plot_isoform_reduced_dim}
\title{FLAMES isoform reduced dimensions plots}
\usage{
plot_isoform_reduced_dim(
  sce,
  gene_id,
  transcript_ids,
  n = 4,
  reduced_dim_name = "UMAP",
  use_gene_dimred = FALSE,
  expr_func = function(x) {
     SingleCellExperiment::logcounts(x)
 },
  col_low = "#313695",
  col_mid = "#FFFFBF",
  col_high = "#A50026",
  color_quantile = 1,
  format = "plot_grid",
  ...
)
}
\arguments{
\item{sce}{The \code{SingleCellExperiment} object containing transcript counts,
\code{rowRanges} and \code{rowData} with \code{gene_id} and \code{transcript_id} columns.}

\item{gene_id}{The gene symbol of interest, ignored if \code{transcript_ids} is provided.}

\item{transcript_ids}{The transcript ids to plot.}

\item{n}{The number of top isoforms to plot from the gene. Ignored if \code{transcript_ids} is provided.}

\item{reduced_dim_name}{The name of the reduced dimension to use for plotting cells.}

\item{use_gene_dimred}{Whether to use gene-level reduced dimensions for plotting. Set to \code{TRUE}
if the \code{SingleCellExperiment} has gene counts in main assay and transcript counts in \code{altExp}.}

\item{expr_func}{The function to extract expression values from the \code{SingleCellExperiment} object.
Default is \code{logcounts}. Alternatively, \code{counts} can be used for raw counts.}

\item{col_low}{Color for cells with low expression levels in UMAPs.}

\item{col_mid}{Color for cells with intermediate expression levels in UMAPs.}

\item{col_high}{Color for cells with high expression levels in UMAPs.}

\item{color_quantile}{The lower and upper expression quantile to be displayed bewteen \code{col_low} and \code{col_high}, e.g. with \code{color_quantile = 0.95}, cells with expressions higher than 95\% of other cells will all be shown in \code{col_high}, and cells with expression lower than 95\% of other cells will all be shown in \code{col_low}.}

\item{format}{The format of the output, either "plot_grid" or "list".}

\item{...}{Additional arguments to pass to \code{plot_grid}.}
}
\value{
a \code{ggplot} object of the UMAP(s)
}
\description{
Plot expression of top n isoforms of a gene in reduced dimensions
}
\details{
Takes \code{SingleCellExperiment} object and plots an expression on reduced dimensions
with the isoform visualizations along genomic coordinates.
}
\examples{
data(scmixology_lib10_transcripts, scmixology_lib10, scmixology_lib90)
scmixology_lib10 <- 
  scmixology_lib10[, colSums(SingleCellExperiment::counts(scmixology_lib10)) > 0]
sce_lr <- scmixology_lib10[, colnames(scmixology_lib10) \%in\% colnames(scmixology_lib10_transcripts)]
SingleCellExperiment::altExp(sce_lr, "transcript") <-
  scmixology_lib10_transcripts[, colnames(sce_lr)]
combined_sce <- combine_sce(sce_lr, scmixology_lib90)
combined_sce <- combined_sce |>
  scuttle::logNormCounts() |>
  scater::runPCA() |>
  scater::runUMAP()
combined_imputed_sce <- sc_impute_transcript(combined_sce)
plot_isoform_reduced_dim(combined_sce, 'ENSG00000108107')
plot_isoform_reduced_dim(combined_imputed_sce, 'ENSG00000108107')

}
