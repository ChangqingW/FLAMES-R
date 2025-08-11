#' FLAMES Differential Transcript Usage Analysis
#'
#' Differential transcription usage testing for single cell data, using
#' \code{colLabels} as cluster labels.
#'
#' @details
#' Genes with more than 2 isoforms expressing more than \code{min_count} counts
#' are selected for testing with one of the following methods:
#' \describe{
#'  \item{trascript usage permutation}{ Transcript usage are taken as the test statistic, cluster labels are permuted to generate a null distribution.}
#'  \item{chisq}{ Chi-square test of the transcript count matrix for each gene. }
#' }
#' Adjusted P-values were calculated by Benjamini–Hochberg correction.
#'
#' @param sce The \code{SingleCellExperiment} object, with transcript counts
#' in the \code{counts} slot and cluster labels in the \code{colLabels} slot.
#' @param min_count The minimum total counts for a transcript to be tested.
#' @param gene_col The column name in the \code{rowData} slot of \code{sce}
#' that contains the gene ID / name. Default is \code{"gene_id"}.
#' @param method The method to use for testing, listed in \code{details}.
#' @param permuations Number of permutations for permutation methods.
#' @param threads Number of threads to use for parallel processing.
#'
#' @return a \code{tibble} containing the following columns:
#' \describe{
#'  \item{p.value}{ - the raw p-value }
#'  \item{adj.p.value}{ - multiple testing adjusted p-value }
#'  \item{cluster}{ - the cluster where DTU was observed }
#'  \item{transcript}{ - rowname of \code{sce}, the DTU isoform }
#'  \item{transcript_usage}{ - the transcript usage of the isoform in the cluster }
#' }
#' Additional columns from \code{method = "trascript usage permutation"}:
#' \describe{
#'  \item{transcript_usage_elsewhere}{ - transcript usage in other clusters }
#'  \item{usage_difference}{ - the difference between the two transcript usage }
#'  \item{permuted_var}{ - the variance of usage difference in the permuted data }
#' }
#' Additional columns from \code{method = "chisq"}:
#' \describe{
#'  \item{X_value}{ - the test statistic }
#'  \item{df}{ - the degrees of freedom }
#'  \item{expected_usage}{ - the expected usage (mean across all clusters) }
#'  \item{usage_difference}{ - the difference between the observed and expected usage }
#' }
#' The table is sorted by P-values.
#'
#' @importFrom methods is
#' @importFrom SingleCellExperiment colLabels
#' @export
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' bc_allow <- file.path(outdir, "bc_allow.tsv")
#' genome_fa <- file.path(outdir, "rps24.fa")
#' R.utils::gunzip(
#'   filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
#'   destname = bc_allow, remove = FALSE
#' )
#' R.utils::gunzip(
#'   filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
#'   destname = genome_fa, remove = FALSE
#' )
#'
#' sce <- FLAMES::sc_long_pipeline(
#'   genome_fa = genome_fa,
#'   fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
#'   annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
#'   outdir = outdir,
#'   barcodes_file = bc_allow,
#'   config_file = create_config(
#'     outdir,
#'     pipeline_parameters.demultiplexer = "flexiplex"
#'   )
#' )
#' group_anno <- data.frame(barcode_seq = colnames(sce), groups = SingleCellExperiment::counts(sce)["ENSMUST00000169826.2", ] > 1)
#' SingleCellExperiment::colLabels(sce) <- group_anno$groups
#' # DTU with permutation testing:
#' sc_DTU_analysis(sce, min_count = 1, method = "trascript usage permutation")
#' # now try with chisq:
#' sc_DTU_analysis(sce, min_count = 1, method = "chisq")
sc_DTU_analysis <- function(sce, gene_col = "gene_id", min_count = 15, threads = 1, method = "trascript usage permutation", permuations = 1000) {

  stopifnot("sce need to be an SingleCellExperiment Object" = is(sce, "SingleCellExperiment"))
  stopifnot("min_count must be a positive number" = min_count > 0)
  stopifnot("Cluster label (colLabels(sce)) not found" = !is.null(colLabels(sce)))
  if (method == "trascript usage permutation") {
    return(
      sc_transcript_usage_permutation(
        sce = sce,
        gene_col = gene_col,
        min_count = min_count,
        threads = threads,
        permuations = permuations
      )
    )
  } else if (method == "chisq") {
    return(
      sc_transcript_usage_chisq(
        sce = sce,
        gene_col = gene_col,
        min_count = min_count,
        threads = threads
      )
    )
  } else {
    stop("Unknown method. Currently available methods are: 'trascript usage permutation' and 'chisq'")
  }
}

#' @importFrom SingleCellExperiment counts colLabels
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr summarise filter n bind_rows arrange mutate
#' @importFrom tibble as_tibble
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom S4Arrays colsum rowsum
#' @keywords internal
sc_transcript_usage_chisq <- function(sce, gene_col = "gene_id", min_count = 15, threads = 1) {

  message("Filtering for genes with at least 2 detected isforms ...")
  sce <- sce[rowSums(SingleCellExperiment::counts(sce)) > min_count, ]
  # TODO: preserve correct gene counts after filtering
  # like in sc_transcript_usage_permutation
  genes <- SummarizedExperiment::rowData(sce) |>
    as.data.frame() |>
    tibble::as_tibble(rownames = "transcript") |>
    dplyr::mutate(n = dplyr::n(), .by = eval(gene_col)) |>
    dplyr::filter(n > 1)
  sce <- sce[genes$transcript, ]
  message(sprintf("\t%d isoform(s) left.\n", nrow(sce)))

  message("Aggregating counts by cluster labels ...")
  pseudobulk_counts <- SingleCellExperiment::counts(sce) |>
    as("CsparseMatrix") |>
    S4Arrays::colsum(
      group = SingleCellExperiment::colLabels(sce)
    )
  split(1:nrow(sce), SummarizedExperiment::rowData(sce)[, gene_col]) |>
    BiocParallel::bplapply(
      FUN = \(x) chisq_test_by_gene(pseudobulk_counts[x, ]),
      BPPARAM = BiocParallel::MulticoreParam(
        workers = threads, stop.on.error = TRUE, progressbar = TRUE
      )
    ) |>
    dplyr::bind_rows(.id = "gene") |>
    dplyr::mutate(adj.p.value = p.adjust(p.value, method = "BH")) |>
    dplyr::arrange(adj.p.value)
}

# a chisq.test wrapper for a gene matrix
chisq_test_by_gene <- function(gene_mtx) {

  # warned <- FALSE
  fit <- withCallingHandlers(
    chisq.test(gene_mtx),
    warning = function(w) {
      if (w$message == "Chi-squared approximation may be incorrect") {
        # warned <<- TRUE
        invokeRestart("muffleWarning")  # Suppress the warning output
      }
      # don't muffle other warnings
  })

  DTU_idx <- arrayInd(which.max(abs(fit$residuals^2)), dim(gene_mtx))
  observed_pcts <- sweep(gene_mtx, 2, colSums(gene_mtx), "/")
  expected_pcts <- rowSums(gene_mtx) / sum(gene_mtx)
  pct_diff <- observed_pcts[DTU_idx[1], DTU_idx[2]] - expected_pcts[DTU_idx[1]]

  tibble::tibble(
    X_value = fit$statistic,
    df = fit$parameter,
    transcript = rownames(gene_mtx)[DTU_idx[1]],
    cluster = colnames(gene_mtx)[DTU_idx[2]],
    p.value = fit$p.value,
    expected_usage = expected_pcts[DTU_idx[1]],
    transcript_usage = observed_pcts[DTU_idx[1], DTU_idx[2]],
    usage_difference = pct_diff
  )
}


# DTU via mean difference with permutation
#' @importFrom MatrixGenerics rowSums
#' @importFrom SingleCellExperiment counts rowData colLabels
#' @importFrom abind abind
#' @importFrom dplyr summarise filter n bind_rows left_join mutate distinct
#' @importFrom BiocParallel bpmapply MulticoreParam
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyselect matches
#' @importFrom S4Arrays rowsum colsum
#' @importFrom stats p.adjust
#' @keywords internal
sc_transcript_usage_permutation <- function(sce, gene_col = "gene_id", min_count = 50, threads = 1, permuations = 1000) {
  message(
    sprintf(
      "Filtering for genes with at least 2 isforms expressing more than %d counts ...",
      min_count
    )
  )
  SummarizedExperiment::rowData(sce)$total <- rowSums(SingleCellExperiment::counts(sce))
  # metadata tibble
  # keep genes with at least 2 isoforms expressing more than min_count
  # all isoforms of a kept gene are kept to get the correct transcript usage (total gene expression)
  genes <- SummarizedExperiment::rowData(sce) |>
    as.data.frame() |>
    dplyr::mutate(test = total >= min_count) |>
    dplyr::mutate(n = sum(test), .by = eval(gene_col)) |>
    dplyr::filter(n > 1) |>
    tibble::as_tibble(rownames = "transcript")
  # filter the sce to remove filtered genes
  sce <- sce[genes$transcript, ]
  message(
    sprintf(
      "\t%d gene(s), %d transcript(s) left.",
      nrow(dplyr::distinct(genes[, gene_col])),
      sum(genes$test)
    )
  )

  cell_labels <- SingleCellExperiment::colLabels(sce)
  gene_ids <- SummarizedExperiment::rowData(sce)[, gene_col]
  gene_counts <- SingleCellExperiment::counts(sce) |>
    as("CsparseMatrix") |> # fixed in SparseArray 1.7.7
    S4Arrays::rowsum(group = gene_ids, reorder = TRUE)
  # with gene_counts saved, can remove filtered isoforms now
  genes <- dplyr::filter(genes, test)
  sce <- sce[genes$transcript, ]
  # update the gene_ids since the filtered isoforms are removed
  gene_ids <- SummarizedExperiment::rowData(sce)[, gene_col]

  transcript_list <- split(rownames(sce), gene_ids)
  # reduce multi-threading overhead by partitioning the job
  # into 4 job partitions per thread
  job_list <- split(transcript_list, cut(seq_along(transcript_list), threads * 4))
  job_list <- job_list[lapply(job_list, length) > 0] # incase of empty job from over partitioning
  # need the gene names since each partition may have different genes
  job_gene_list <- lapply(job_list, function(jbs) {
    lapply(names(jbs), \(gene) rep(gene, length(jbs[[gene]]))) |>
      do.call(what = c, args = _)
  })
  job_list <- lapply(job_list, unlist)

  tb <- BiocParallel::bpmapply(
    function(transcripts, genes) {
      mtx <- SingleCellExperiment::counts(sce)[transcripts, ] |>
        as("CsparseMatrix")
      orig_diff_mtx <- mean_transcript_usage(mtx, cell_labels, genes, gene_counts[genes, ])

      # permute the labels and get the mean difference matrix
      perm_diff_mtx <- lapply(seq_len(permuations), \(x) {
        perm_labels <- sample(cell_labels, length(cell_labels), replace = FALSE)
        mean_transcript_usage(mtx, perm_labels, genes, gene_counts[genes, ], diff_only = TRUE)
      }) |>
        abind::abind(along = 3) |>
        abs()


      args.grid <- expand.grid(
        transcript = rownames(orig_diff_mtx),
        cluster = colnames(orig_diff_mtx)
      )
      mapply(
        function(i, j) {
          permuted_diffs <- perm_diff_mtx[i, j, ]
          n_samples <- length(na.omit(permuted_diffs))
          # in case of ties, get the lower bound of the quantile
          imperical_quantile <- ecdf_lower(permuted_diffs, abs(orig_diff_mtx[i, j, "usage_difference"]))
          p_value <- (1 - imperical_quantile) + (1 / n_samples) # percision at best is 1 /n_samples 
          permuted_var <- sum(permuted_diffs^2, na.rm = TRUE) / (sum(!is.na(permuted_diffs)) - 1)
          return(tibble::tibble(
            transcript = i,
            cluster = j,
            transcript_usage = orig_diff_mtx[i, j, "transcript_usage"],
            transcript_usage_elsewhere = orig_diff_mtx[i, j, "transcript_usage_elsewhere"],
            usage_difference = orig_diff_mtx[i, j, "usage_difference"],
            p.value = p_value,
            # permuted_mean = permuted_mean,
            permuted_var = permuted_var
          ))
        },
        args.grid$transcript, args.grid$cluster,
        SIMPLIFY = FALSE
      ) |>
        dplyr::bind_rows()
    },
    job_list, job_gene_list,
    SIMPLIFY = FALSE,
    BPPARAM = BiocParallel::MulticoreParam(workers = threads, progressbar = TRUE)
  ) |>
    dplyr::bind_rows()

  # sanity check
  stopifnot("Unexpected error" = nrow(tb) == (nrow(sce) * length(unique(cell_labels))))
  rowdata <- SummarizedExperiment::rowData(sce) |>
    as.data.frame() |>
    tibble::as_tibble(rownames = "transcript") |>
    dplyr::select(transcript, eval(gene_col), tidyselect::matches("gene_name"))
  dplyr::mutate(tb, adj.p.value = p.adjust(p.value, method = "BH")) |>
    dplyr::left_join(rowdata, by = c("transcript")) |>
    dplyr::arrange(adj.p.value)
}

ecdf_lower <- function(samples, x) {
  sum(samples < x) / length(samples)
}

#' @importFrom S4Arrays rowsum colsum
#' @importFrom MatrixGenerics rowSums
#' @importFrom abind abind
#' @keywords internal
mean_transcript_usage <- function(mtx, cell_labels, genes, gene_counts, diff_only = FALSE) {
  stopifnot("labels should be of the same length as the number of columns in mtx" = length(cell_labels) == ncol(mtx))

  transcript_counts <- mtx |>
    S4Arrays::colsum(group = cell_labels, reorder = TRUE)
  if (missing(genes)) {
    genes <- rep(1, nrow(mtx))
  } else {
    stopifnot("genes should be of the same length as the number of rows in mtx" = length(genes) == nrow(mtx))
  }

  if (missing(gene_counts)) {
    # if trascript counts are filtered, the transcript usage will be wrong!
    gene_counts <- transcript_counts |>
      S4Arrays::rowsum(group = genes, reorder = TRUE) |>
      (\(mtx) mtx[genes, ])() # need same number of rows as transcript count mtx
  } else {
    stopifnot("gene_counts have the same dimension as mtx" = all(dim(gene_counts) == dim(mtx)))
    gene_counts <- gene_counts |>
      S4Arrays::colsum(group = cell_labels, reorder = TRUE)
  }
  transcript_usage <- transcript_counts / gene_counts

  # get the mean transcript usage of other clusters (ecxluding the current cluster)
  transcript_usage_elsewhere <- (rowSums(transcript_counts) - transcript_counts) / (rowSums(gene_counts) - gene_counts)

  usage_difference <- transcript_usage - transcript_usage_elsewhere
  if (diff_only) {
    return(usage_difference)
  }

  diff_mtx <- abind::abind(transcript_usage,
    transcript_usage_elsewhere,
    usage_difference,
    along = 3
  )
  dimnames(diff_mtx)[[3]] <- c("transcript_usage", "transcript_usage_elsewhere", "usage_difference")
  return(diff_mtx)
}
