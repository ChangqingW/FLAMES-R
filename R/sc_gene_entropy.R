# Modified from Sefi Prawer's scripts
# sefi.prawer@unimelb.edu.au
#' Compute Gene Isoform Entropy Matrix
#'
#' Calculates normalized Shannon entropy for gene isoform expression across cells.
#' Higher entropy indicates more diverse isoform usage, lower entropy indicates
#' dominance by fewer isoforms.
#'
#' @param sce A \code{SingleCellExperiment} object
#' @param assay Name of assay containing isoform counts (default: "counts")
#' @param gene_col Column name in rowData containing gene identifiers (default: "gene_id")
#' @param alpha Pseudocount added to avoid log(0) (default: .Machine$double.xmin)
#' @param min_counts_per_cell Minimum total gene counts per cell to include (default: 5)
#' @param isoform_min_pct_cells Minimum fraction of cells expressing each isoform (default: 0.05)
#' @param isoform_cumulative_pct Keep top isoforms contributing to this cumulative proportion (default: 0.95)
#' @param min_cell_fraction Minimum fraction of cells with valid entropy per gene (default: 0.25)
#' @param threads Number of threads for parallel processing (default: 1)
#' @param show_progress Logical indicating whether to show progress (default: TRUE if interactive)
#'
#' @return Matrix with genes as rows and cells as columns containing normalized entropy values (0-1).
#'
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom cli cli_alert_info
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom Matrix sparseMatrix
#' @examples
#' sce <- scuttle::mockSCE(ncells = 50, ngenes = 30)
#' SummarizedExperiment::rowData(sce)$gene_id <- sort(
#'   paste0("gene", sample(1:9, nrow(sce), replace = TRUE))
#' )
#' res <- sc_gene_entropy(sce, threads = 2)
#'
#' @export
sc_gene_entropy <- function(
    sce, assay = "counts", gene_col = "gene_id", alpha = .Machine$double.xmin,
    min_counts_per_cell = 5, isoform_min_pct_cells = 0.05,
    isoform_cumulative_pct = 0.95, min_cell_fraction = 0.25, threads = 1,
    show_progress = interactive()) {
  if (!methods::is(sce, "SingleCellExperiment")) {
    stop("sce must be a SingleCellExperiment object")
  }

  if (show_progress) {
    cli::cli_alert_info("Pre-processing data and building gene-isoform mapping...")
  }

  cells <- colnames(sce)
  mat_counts <- SummarizedExperiment::assay(sce, assay)[, cells, drop = FALSE]

  # Convert to sparse matrix for memory efficiency
  if (!inherits(mat_counts, "sparseMatrix")) {
    mat_counts <- as(mat_counts, "CsparseMatrix")
  }

  # Pre-build gene-to-isoform mapping to avoid repeated grep operations
  rowdata <- SummarizedExperiment::rowData(sce)
  if (gene_col %in% colnames(rowdata)) {
    # Use existing gene annotation
    gene_to_isoforms <- split(rownames(mat_counts), rowdata[[gene_col]])
  } else {
    stop(paste("Column", gene_col, "not found in rowData. Please provide a valid gene identifier column."))
  }

  # Pre-filter genes with <2 isoforms
  valid_genes <- names(gene_to_isoforms)[lengths(gene_to_isoforms) >= 2]
  gene_to_isoforms <- gene_to_isoforms[valid_genes]

  if (length(gene_to_isoforms) == 0) {
    stop("No genes found with 2 or more isoforms")
  }

  if (show_progress) {
    cli::cli_alert_info("Processing {length(gene_to_isoforms)} genes with {length(cells)} cells...")
  }

  # Reduce multi-threading overhead by partitioning jobs into chunks
  # Use 4 job partitions per thread to balance overhead vs parallelism
  n_jobs <- max(1, threads * 4)
  job_chunks <- split(seq_along(gene_to_isoforms), cut(seq_along(gene_to_isoforms), n_jobs))
  job_chunks <- job_chunks[lengths(job_chunks) > 0]

  # Process chunks in parallel with BiocParallel's progress bar
  chunk_results <- BiocParallel::bplapply(
    job_chunks,
    function(chunk_indices) {
      chunk_genes <- names(gene_to_isoforms)[chunk_indices]
      chunk_mapping <- gene_to_isoforms[chunk_indices]

      # Pre-allocate result matrix for this chunk
      chunk_entropy_mat <- matrix(NA_real_,
        nrow = length(chunk_genes),
        ncol = length(cells)
      )
      rownames(chunk_entropy_mat) <- chunk_genes
      colnames(chunk_entropy_mat) <- cells

      # Process genes in this chunk
      for (i in seq_along(chunk_genes)) {
        gene <- chunk_genes[i]
        isoforms <- chunk_mapping[[i]]

        # Extract and process data for this gene
        submat <- mat_counts[isoforms, , drop = FALSE]
        gene_total_counts <- Matrix::colSums(submat)
        keep_cells <- names(gene_total_counts)[gene_total_counts >= min_counts_per_cell]

        if (length(keep_cells) == 0) next

        submat_f <- submat[, keep_cells, drop = FALSE]

        # Vectorized calculation of expressed isoforms percentage
        pct_cells_expressed <- Matrix::rowMeans(submat_f > 0)
        iso_keep_pct <- names(pct_cells_expressed)[pct_cells_expressed >= isoform_min_pct_cells]
        submat_f <- submat_f[iso_keep_pct, , drop = FALSE]

        if (nrow(submat_f) < 2) next

        # Optimized top isoform filtering using vectorized operations
        col_totals <- Matrix::colSums(submat_f)
        prop_mat <- sweep(as.matrix(submat_f), 2, col_totals, FUN = "/")

        # Find top isoforms per cell more efficiently
        iso_keep_final <- character(0)
        for (j in seq_len(ncol(prop_mat))) {
          if (col_totals[j] == 0) next
          sorted_props <- sort(prop_mat[, j], decreasing = TRUE)
          cum_props <- cumsum(sorted_props)
          n_keep <- max(1, sum(cum_props <= isoform_cumulative_pct) + 1)
          top_isos <- names(sort(prop_mat[, j], decreasing = TRUE)[seq_len(n_keep)])
          iso_keep_final <- union(iso_keep_final, top_isos)
        }

        iso_keep_final <- intersect(rownames(submat_f), iso_keep_final)
        if (length(iso_keep_final) < 2) next

        submat_final <- submat_f[iso_keep_final, , drop = FALSE]

        # Vectorized entropy calculation
        mat_smoothed <- as.matrix(submat_final) + alpha
        col_sums <- colSums(mat_smoothed)
        p_mat <- sweep(mat_smoothed, 2, col_sums, FUN = "/")

        # Efficient entropy calculation
        log_p_mat <- log2(p_mat)
        log_p_mat[!is.finite(log_p_mat)] <- 0 # Handle log(0)

        entropy_vec <- -colSums(p_mat * log_p_mat) / log2(nrow(p_mat))

        chunk_entropy_mat[gene, colnames(submat_final)] <- entropy_vec
      }

      return(chunk_entropy_mat)
    },
    BPPARAM = BiocParallel::MulticoreParam(workers = threads, progressbar = show_progress)
  )

  # Combine results from all chunks
  entropy_mat <- do.call(rbind, chunk_results)

  # Filter genes with sufficient cell coverage
  entropy_mat <- entropy_mat[rowMeans(!is.na(entropy_mat)) >= min_cell_fraction, , drop = FALSE]

  if (show_progress) {
    cli::cli_alert_info("Completed processing. Final matrix: {nrow(entropy_mat)} genes x {ncol(entropy_mat)} cells")
  }

  return(entropy_mat)
}
