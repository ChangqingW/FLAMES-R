make_dtu_sce <- function(seed = 42) {
  set.seed(seed)
  n_cells <- 20
  # 3 genes x 2 transcripts each; ensure counts are high enough for min_count = 1
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

test_that("sc_DTU_analysis chisq method returns a tibble with expected columns", {
  sce <- make_dtu_sce()
  res <- sc_DTU_analysis(sce, min_count = 1, method = "chisq")
  expect_s3_class(res, "tbl_df")
  expect_true(all(c("p.value", "adj.p.value", "transcript", "cluster") %in% names(res)))
})

test_that("sc_DTU_analysis permutation method returns a tibble with expected columns", {
  sce <- make_dtu_sce()
  res <- sc_DTU_analysis(sce, min_count = 1,
                         method = "transcript usage permutation", permuations = 20)
  expect_s3_class(res, "tbl_df")
  expect_true(all(c("p.value", "adj.p.value", "transcript", "cluster",
                    "transcript_usage", "usage_difference") %in% names(res)))
})

test_that("sc_DTU_analysis rejects non-SCE input", {
  expect_error(sc_DTU_analysis(data.frame()), "SingleCellExperiment")
})

test_that("sc_DTU_analysis errors when colLabels are missing", {
  sce <- make_dtu_sce()
  SingleCellExperiment::colLabels(sce) <- NULL
  expect_error(sc_DTU_analysis(sce, min_count = 1), "Cluster label")
})

test_that("sc_DTU_analysis errors on unknown method", {
  sce <- make_dtu_sce()
  expect_error(sc_DTU_analysis(sce, min_count = 1, method = "bad_method"),
               "Unknown method")
})
