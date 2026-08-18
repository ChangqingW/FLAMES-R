# Tests for the pure-R output parsers create_sce_from_dir / create_se_from_dir
# (+ generate_sc_sce / addRowRanges), driven by small fabricated output dirs so
# no pipeline run / external tool is needed. Transcript IDs are real Rps24 IDs so
# addRowRanges matches the bundled annotation.

.rps24_gtf <- function() system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
.rps24_tx <- c("ENSMUST00000169826.2", "ENSMUST00000225994.1", "ENSMUST00000225023.1")

# Write a FLAMES-style output dir: <prefix>transcript_count.csv.gz +
# isoform_FSM_annotation.csv (+ optional gene_count.* for add_gene_counts).
write_flames_outdir <- function(dir, prefix = "", gene_counts = TRUE) {
  counts <- data.frame(
    transcript_id = .rps24_tx,
    gene_id = "ENSMUSG00000025290.17",
    cell1 = c(5, 2, 1), cell2 = c(3, 0, 4),
    stringsAsFactors = FALSE
  )
  readr::write_csv(counts, file.path(dir, paste0(prefix, "transcript_count.csv.gz")))
  fsm <- data.frame(transcript_id = .rps24_tx, FSM_match = .rps24_tx, stringsAsFactors = FALSE)
  readr::write_csv(fsm, file.path(dir, "isoform_FSM_annotation.csv"))
  if (gene_counts) {
    gm <- Matrix::Matrix(c(1, 2, 3, 4), nrow = 2, sparse = TRUE,
                         dimnames = list(c("geneA", "geneB"), c("cell1", "cell2")))
    Matrix::writeMM(gm, file.path(dir, "gene_count.mtx"))
    writeLines(rownames(gm), file.path(dir, "gene_count_features.tsv"))
    writeLines(colnames(gm), file.path(dir, "gene_count_barcodes.tsv"))
  }
  dir
}

# Write an Oarfish single-cell output dir: <stem>.count.mtx/.features.txt/.barcodes.txt
write_oarfish_sc_outdir <- function(dir, stem = "sample") {
  m <- Matrix::Matrix(c(3, 0, 5, 2, 0, 7), nrow = 2, ncol = 3, sparse = TRUE)
  Matrix::writeMM(m, file.path(dir, paste0(stem, ".count.mtx")))
  writeLines(.rps24_tx, file.path(dir, paste0(stem, ".features.txt")))
  writeLines(c("AAACCCAAGAAACCCA", "AAACCCAAGAAACGGT"), file.path(dir, paste0(stem, ".barcodes.txt")))
  dir
}

# Write an Oarfish bulk output dir: <sample>.quant TSVs
write_oarfish_bulk_outdir <- function(dir) {
  for (s in c("s1", "s2")) {
    writeLines(c("tname\tlen\tnum_reads",
                 paste(.rps24_tx, c(1000, 2000, 1500), c(5, 10, 0), sep = "\t")),
               file.path(dir, paste0(s, ".quant")))
  }
  dir
}

# ---- create_sce_from_dir ---------------------------------------------------

test_that("create_sce_from_dir parses a FLAMES output dir into an SCE with gene counts", {
  dir <- write_flames_outdir(withr::local_tempdir())
  sce <- create_sce_from_dir(dir, .rps24_gtf())
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(ncol(sce), 2L)
  expect_true("gene" %in% SingleCellExperiment::altExpNames(sce))
  expect_equal(length(SummarizedExperiment::rowRanges(sce)), nrow(sce))
})

test_that("create_sce_from_dir parses an Oarfish output dir", {
  dir <- write_oarfish_sc_outdir(withr::local_tempdir())
  # missing quantification -> auto-selects Oarfish since no FLAMES counts present
  expect_message(
    sce <- create_sce_from_dir(dir, .rps24_gtf()),
    "No novel annotation found|Gene counts not added"
  )
  expect_s4_class(sce, "SingleCellExperiment")
})

test_that("create_sce_from_dir errors when the folder has no results", {
  empty <- withr::local_tempdir()
  expect_error(create_sce_from_dir(empty, .rps24_gtf()), "No transcript count results")
  expect_error(create_sce_from_dir(empty, .rps24_gtf(), quantification = "FLAMES"),
               "transcript_count.csv.gz")
  expect_error(create_sce_from_dir(empty, .rps24_gtf(), quantification = "Oarfish"),
               "count.mtx")
})

# ---- create_se_from_dir ----------------------------------------------------

test_that("create_se_from_dir parses a FLAMES output dir into a SummarizedExperiment", {
  dir <- write_flames_outdir(withr::local_tempdir(), gene_counts = FALSE)
  se <- create_se_from_dir(dir, .rps24_gtf())
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(ncol(se), 2L)
})

test_that("create_se_from_dir parses an Oarfish output dir (auto-detected)", {
  dir <- write_oarfish_bulk_outdir(withr::local_tempdir())
  se <- create_se_from_dir(dir, .rps24_gtf())
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(ncol(se), 2L)
})

test_that("create_se_from_dir errors on an unknown quantification method", {
  dir <- write_flames_outdir(withr::local_tempdir(), gene_counts = FALSE)
  expect_error(create_se_from_dir(dir, .rps24_gtf(), quantification = "BOGUS"),
               "Unknown quantification method")
})

# ---- addRowRanges direct branches ------------------------------------------

.mk_sce <- function(rn) {
  m <- matrix(1L, nrow = length(rn), ncol = 2, dimnames = list(rn, c("c1", "c2")))
  SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))
}

test_that("addRowRanges warns for transcripts absent from the annotation", {
  dir <- withr::local_tempdir()
  sce <- .mk_sce(c(.rps24_tx[1], "FAKE_TX"))
  expect_warning(
    res <- FLAMES:::addRowRanges(sce, .rps24_gtf(), dir),
    "not recorded in the annotation"
  )
  expect_s4_class(res, "SingleCellExperiment")
})

test_that("addRowRanges merges a novel isoform_annotated.gtf from the outdir", {
  dir <- withr::local_tempdir()
  R.utils::gunzip(.rps24_gtf(), destname = file.path(dir, "isoform_annotated.gtf"),
                  remove = FALSE, overwrite = TRUE)
  sce <- .mk_sce(.rps24_tx)
  res <- FLAMES:::addRowRanges(sce, .rps24_gtf(), dir)
  expect_equal(length(SummarizedExperiment::rowRanges(res)), nrow(res))
})
