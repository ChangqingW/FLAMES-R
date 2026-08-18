test_that("get_coverage returns per-transcript binned coverage from a BAM", {
  bam <- local_bam()
  cov <- get_coverage(bam, min_counts = 1)
  expect_s3_class(cov, "tbl_df")
  expect_true(all(c("transcript", "read_counts", "tr_length") %in% names(cov)))
  expect_true(any(grepl("^coverage_", names(cov))))
  expect_gte(nrow(cov), 1)
  # chr1 carried 4 reads (>= min_counts)
  expect_true("chr1" %in% cov$transcript)
})

test_that("plot_coverage returns a ggplot for a coverage table", {
  bam <- local_bam()
  cov <- get_coverage(bam, min_counts = 1)
  expect_s3_class(suppressWarnings(plot_coverage(cov)), "ggplot")
})

test_that("filter_annotation returns a GRanges and reduces transcript count", {
  gtf <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
  for (keep_val in c("tss_differ", "tes_differ", "both")) {
    res <- filter_annotation(gtf, keep = keep_val)
    expect_s4_class(res, "GRanges")
    full <- rtracklayer::import(gtf, feature.type = "transcript")
    expect_lte(length(res), length(full))
  }
})

test_that("filter_annotation keep = 'single_transcripts' returns only sole-transcript genes", {
  gtf <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
  res <- filter_annotation(gtf, keep = "single_transcripts")
  expect_s4_class(res, "GRanges")
  # every gene_id in the result must appear exactly once
  gene_counts <- table(S4Vectors::mcols(res)$gene_id)
  expect_true(all(gene_counts == 1))
})

test_that("filter_annotation errors on unknown keep value", {
  gtf <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
  expect_error(filter_annotation(gtf, keep = "invalid"), "Unknown keep value")
})

test_that("weight_transcripts type = 'counts' returns input unchanged", {
  x <- c(5, 50, 500, 5000)
  expect_equal(weight_transcripts(x, type = "counts"), x)
})

test_that("weight_transcripts type = 'equal' produces binary weights", {
  w <- weight_transcripts(c(500, 1001, 999, 2000), type = "equal", min_counts = 1000)
  expect_equal(w, c(0, 1, 0, 1))
})

test_that("weight_transcripts type = 'sigmoid' returns values in [0, 1]", {
  w <- weight_transcripts(1:200, type = "sigmoid")
  expect_true(all(w >= 0 & w <= 1))
  # weights should increase with count
  expect_true(all(diff(w) >= 0))
})

test_that("convolution_filter passes flat coverage and fails sharp drops", {
  # flat coverage -> no abrupt change -> should pass
  expect_true(convolution_filter(rep(1.0, 100)))
  # 50% drop in the middle -> should fail
  expect_false(convolution_filter(c(rep(1.0, 50), rep(0.5, 50))))
  # examples from docstring
  expect_false(convolution_filter(c(1, 1, 1, 0.69, 0.69, 0.69), threshold = 0.3))
  expect_true(convolution_filter(c(1, 1, 1, 0.71, 0.7, 0.7), threshold = 0.3))
})

test_that("filter_coverage retains passing transcripts and removes failing ones", {
  # Build a minimal coverage tibble
  n_pos <- 100
  cov_df <- data.frame(
    transcript = c("pass_t", "fail_t"),
    read_counts = c(100L, 100L),
    tr_length   = c(1000L, 1000L),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(n_pos)) {
    cov_df[[paste0("coverage_", i)]] <- c(
      1.0,                                     # pass_t: flat
      if (i <= 50) 1.0 else 0.5               # fail_t: sharp 50% drop
    )
  }

  result <- filter_coverage(cov_df)
  expect_equal(nrow(result), 1L)
  expect_equal(result$transcript, "pass_t")
})

test_that("plot_coverage(detailed = TRUE) returns a composite plot", {
  cov <- get_coverage(local_bam(), min_counts = 1)
  p <- suppressWarnings(plot_coverage(cov, detailed = TRUE))
  expect_true(inherits(p, c("ggplot", "gg", "gtable"))) # cowplot::plot_grid result
})

test_that("plot_coverage supports fixed length bins", {
  cov <- get_coverage(local_bam(), min_counts = 1)
  p <- suppressWarnings(plot_coverage(cov, length_bins = c(0, 1, 2, 5, 10, Inf)))
  expect_s3_class(p, "ggplot")
})

test_that("plot_coverage and filter_coverage accept a BAM path", {
  bam <- local_coverage_bam()
  expect_s3_class(suppressWarnings(plot_coverage(bam, filter_fn = convolution_filter)), "ggplot")
  expect_true(is.data.frame(filter_coverage(bam)))
})

test_that("weight_transcripts sigmoid caps the inflection index for large inputs", {
  w <- weight_transcripts(1:2000, type = "sigmoid")
  expect_length(w, 2000)
  expect_true(all(w >= 0 & w <= 1))
})
