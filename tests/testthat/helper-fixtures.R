# Shared helpers for the FLAMES tests

# barcode allow-list
local_bc_allow <- function(envir = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = envir)
  bc_allow <- file.path(dir, "bc_allow.tsv")
  R.utils::gunzip(
    filename = system.file("extdata/bc_allow.tsv.gz", package = "FLAMES"),
    destname = bc_allow, remove = FALSE, overwrite = TRUE
  )
  bc_allow
}

# The 10x barcode pattern
barcode_pattern <- function() {
  c(
    primer = "CTACACGACGCTCTTCCGATCT",
    BC = paste0(rep("N", 16), collapse = ""),
    UMI = paste0(rep("N", 12), collapse = ""),
    polyT = paste0(rep("T", 9), collapse = "")
  )
}

# A small SingleCellExperiment for DTU testing: 3 genes x 2 transcripts each,
# `n_cells` cells split evenly into two clusters (A / B). Counts are Poisson
# with a high enough mean to clear min_count = 1.
make_dtu_sce <- function(seed = 42, n_cells = 20) {
  withr::local_seed(seed)
  counts_mat <- matrix(
    rpois(6 * n_cells, lambda = 15),
    nrow = 6, ncol = n_cells,
    dimnames = list(paste0("tx", 1:6), paste0("cell", 1:n_cells))
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts_mat)
  )
  SummarizedExperiment::rowData(sce)$gene_id <- rep(paste0("gene", 1:3), each = 2)
  SingleCellExperiment::colLabels(sce) <- factor(rep(c("A", "B"), n_cells / 2))
  sce
}

# A mock SingleCellExperiment carrying a rowData$gene_id column that groups the
# `ngenes` rows into `ngroups` genes, so most genes have >= 2 isoforms (as
# find_diversity requires).
make_diversity_sce <- function(ncells = 50, ngenes = 30, ngroups = 9, seed = 42) {
  withr::local_seed(seed)
  sce <- scuttle::mockSCE(ncells = ncells, ngenes = ngenes)
  SummarizedExperiment::rowData(sce)$gene_id <- sort(
    paste0("gene", sample(seq_len(ngroups), nrow(sce), replace = TRUE))
  )
  sce
}
