test_that("find_diversity returns a numeric matrix with correct shape", {
  set.seed(42)
  sce <- scuttle::mockSCE(ncells = 50, ngenes = 30)
  SummarizedExperiment::rowData(sce)$gene_id <- sort(
    paste0("gene", sample(1:9, nrow(sce), replace = TRUE))
  )

  res <- find_diversity(sce, threads = 1, show_progress = FALSE)

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_equal(ncol(res), ncol(sce))
  expect_true(nrow(res) > 0)
  # values should be NA or in [0, 1]
  expect_true(all(is.na(res) | (res >= 0 & res <= 1)))
})

test_that("find_diversity rejects non-SCE input", {
  expect_error(find_diversity(list()), "SingleCellExperiment")
  expect_error(find_diversity(data.frame()), "SingleCellExperiment")
})

test_that("find_diversity errors when gene_col is absent", {
  sce <- scuttle::mockSCE(ncells = 10, ngenes = 10)
  expect_error(find_diversity(sce, gene_col = "no_such_col", show_progress = FALSE),
               "not found in rowData")
})

test_that("find_diversity errors when no gene has >= 2 isoforms", {
  set.seed(1)
  sce <- scuttle::mockSCE(ncells = 10, ngenes = 4)
  # every row is its own gene -> no gene has 2 isoforms
  SummarizedExperiment::rowData(sce)$gene_id <- paste0("g", seq_len(nrow(sce)))
  expect_error(find_diversity(sce, show_progress = FALSE), "2 or more isoforms")
})
