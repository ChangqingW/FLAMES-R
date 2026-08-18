# test excludes the basilisk/python, minimap2, or oarfish side

test_that("add_gene_counts adds a gene altExp aligned to the SCE columns", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0, nrow = 10, ncol = 10))
  )
  colnames(sce) <- paste0("cell", 1:10)

  stem <- file.path(withr::local_tempdir(), "gene_count")
  gene_mtx <- Matrix::Matrix(1:10, nrow = 2, ncol = 5, sparse = TRUE)
  colnames(gene_mtx) <- paste0("cell", 1:5)
  rownames(gene_mtx) <- c("gene1", "gene2")
  Matrix::writeMM(gene_mtx, paste0(stem, ".mtx"))
  writeLines(rownames(gene_mtx), paste0(stem, "_features.tsv"))
  writeLines(colnames(gene_mtx), paste0(stem, "_barcodes.tsv"))

  res <- add_gene_counts(sce, stem)
  gene <- SingleCellExperiment::altExps(res)$gene

  expect_s4_class(gene, "SingleCellExperiment")
  expect_equal(dim(gene), c(2, 10))
  expect_identical(colnames(gene), colnames(sce))
  expect_identical(rownames(gene), c("gene1", "gene2"))

  cnts <- as.matrix(SingleCellExperiment::counts(gene))
  # cells 6-10 were absent from the gene matrix, should be zero-filled
  expect_true(all(cnts[, 6:10] == 0))
  # total counts are preserved for the cells that were present
  expect_equal(sum(cnts), sum(1:10))
})

test_that("add_gene_counts errors when the gene count files are missing", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0, nrow = 2, ncol = 2))
  )
  colnames(sce) <- c("cell1", "cell2")
  expect_error(
    add_gene_counts(sce, file.path(tempdir(), "definitely_not_here")),
    "not found"
  )
})

test_that("parse_oarfish_bulk_output builds a SummarizedExperiment of num_reads", {
  dir <- withr::local_tempdir()
  make_quant <- function(stem, num_reads) {
    writeLines(
      c(
        "tname\tlen\tnum_reads",
        paste("tx1", 1000, num_reads[1], sep = "\t"),
        paste("tx2", 2000, num_reads[2], sep = "\t"),
        paste("tx3", 1500, num_reads[3], sep = "\t")
      ),
      paste0(stem, ".quant")
    )
    stem
  }
  s1 <- make_quant(file.path(dir, "s1"), c(5, 10, 0))
  s2 <- make_quant(file.path(dir, "s2"), c(1, 2, 3))

  se <- FLAMES:::parse_oarfish_bulk_output(c(s1, s2), c("sample1", "sample2"))

  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(dim(se), c(3, 2))
  expect_identical(colnames(se), c("sample1", "sample2"))
  expect_identical(rownames(se), c("tx1", "tx2", "tx3"))

  m <- as.matrix(SummarizedExperiment::assay(se, "counts"))
  expect_equal(unname(m[, "sample1"]), c(5, 10, 0))
  expect_equal(unname(m[, "sample2"]), c(1, 2, 3))
})

test_that("parse_oarfish_sc_output builds an SCE with populated rowRanges", {
  dir <- withr::local_tempdir()
  stem <- file.path(dir, "oarfish_sc")
  # transcript IDs from the Rps24 annotation
  features <- c("ENSMUST00000169826.2", "ENSMUST00000225994.1", "ENSMUST00000225023.1")
  barcodes <- c("AAACCCAAGAAACCCA", "AAACCCAAGAAACGGT")
  # oarfish writes a (barcodes x features) matrix; the parser transposes it
  m <- Matrix::Matrix(c(3, 0, 5, 2, 0, 7), nrow = 2, ncol = 3, sparse = TRUE)
  Matrix::writeMM(m, paste0(stem, ".count.mtx"))
  writeLines(features, paste0(stem, ".features.txt"))
  writeLines(barcodes, paste0(stem, ".barcodes.txt"))

  gtf <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
  # no isoform_annotated.* in outdir -> "No novel annotation found."
  expect_message(
    sce <- FLAMES:::parse_oarfish_sc_output(stem, gtf, outdir = dir),
    "No novel annotation found"
  )

  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(dim(sce), c(3, 2))
  expect_identical(rownames(sce), features)
  expect_identical(colnames(sce), barcodes)
  expect_true("transcript_id" %in% colnames(SummarizedExperiment::rowData(sce)))
  expect_equal(length(SummarizedExperiment::rowRanges(sce)), 3L)
})

test_that("quantify_transcript_oarfish errors on an unknown pipeline", {
  gtf <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
  expect_error(
    FLAMES:::quantify_transcript_oarfish(
      gtf, withr::local_tempdir(), config = list(), pipeline = "nope", samples = "s"
    ),
    "Unknown pipeline"
  )
})

test_that("add_gene_counts derives the count stem from metadata outdir", {
  dir <- withr::local_tempdir()
  gm <- Matrix::Matrix(1:10, nrow = 2, ncol = 5, sparse = TRUE)
  colnames(gm) <- paste0("cell", 1:5)
  rownames(gm) <- c("g1", "g2")
  Matrix::writeMM(gm, file.path(dir, "gene_count.mtx"))
  writeLines(rownames(gm), file.path(dir, "gene_count_features.tsv"))
  writeLines(colnames(gm), file.path(dir, "gene_count_barcodes.tsv"))

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0, nrow = 2, ncol = 10))
  )
  colnames(sce) <- paste0("cell", 1:10)
  S4Vectors::metadata(sce)$inputs$outdir <- dir

  res <- add_gene_counts(sce)
  expect_true("gene" %in% SingleCellExperiment::altExpNames(res))
})

test_that("parse_oarfish_sc_output warns for transcripts absent from the annotation", {
  dir <- withr::local_tempdir()
  stem <- file.path(dir, "oarfish")
  m <- Matrix::Matrix(c(3, 0, 5, 2), nrow = 2, ncol = 2, sparse = TRUE)
  Matrix::writeMM(m, paste0(stem, ".count.mtx"))
  writeLines(c("ENSMUST00000169826.2", "FAKE_TX"), paste0(stem, ".features.txt"))
  writeLines(c("bc1", "bc2"), paste0(stem, ".barcodes.txt"))
  gtf <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")

  expect_warning(
    suppressMessages(FLAMES:::parse_oarfish_sc_output(stem, gtf, dir)),
    "not found in the annotation"
  )
})
