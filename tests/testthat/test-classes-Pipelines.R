# Tests for the pipeline classes' show()/display helpers (classes-Pipelines.R).

# helpers

test_that("truncate_path shortens long paths and leaves short ones intact", {
  expect_identical(FLAMES:::truncate_path("/short/path"), "/short/path")
  long <- paste0("/a/very/long/", paste(rep("dir", 40), collapse = "/"), "/file.bam")
  out <- FLAMES:::truncate_path(long, max_chars = 50)
  expect_true(startsWith(out, "..."))
  expect_equal(nchar(out), 50)
})

test_that("format_file_size renders bytes in human units", {
  expect_equal(FLAMES:::format_file_size(0), "0 B")
  expect_equal(FLAMES:::format_file_size(NA), "0 B")
  expect_equal(FLAMES:::format_file_size(512), "512.0 B")
  expect_equal(FLAMES:::format_file_size(1024), "1.0 KB")
  expect_equal(FLAMES:::format_file_size(1024^2), "1.0 MB")
  expect_equal(FLAMES:::format_file_size(3 * 1024^3), "3.0 GB")
})

test_that("get_expect_cell_display handles unset / single / multi / index", {
  ppl <- example_pipeline("SingleCellPipeline")
  ppl@expect_cell_number <- numeric(0)
  expect_null(FLAMES:::get_expect_cell_display(ppl))
  ppl@expect_cell_number <- 100
  expect_equal(FLAMES:::get_expect_cell_display(ppl), 100)
  ppl@expect_cell_number <- c(100, 200, 300)
  expect_equal(FLAMES:::get_expect_cell_display(ppl), "100, 200, 300")
  expect_equal(FLAMES:::get_expect_cell_display(ppl, index = 2), 200)
  expect_null(FLAMES:::get_expect_cell_display(ppl, index = 9))
})

test_that("format_durations picks appropriate units", {
  d <- as.difftime(c(5, 120, 7200), units = "secs")
  out <- FLAMES:::format_durations(d)
  expect_match(out[1], "sec")
  expect_match(out[2], "min")
  expect_match(out[3], "hr")
})

# show()

test_that("show() runs for all three pipeline types (fresh objects)", {
  for (type in c("BulkPipeline", "SingleCellPipeline", "MultiSampleSCPipeline")) {
    ppl <- example_pipeline(type)
    expect_no_error(show(ppl))
  }
})

test_that("show() renders completed, failed and pending steps", {
  ppl <- example_pipeline("SingleCellPipeline")
  ppl@completed_steps[c("barcode_demultiplex", "genome_alignment")] <- TRUE
  d <- as.difftime(c(5, 130), units = "secs")
  names(d) <- c("barcode_demultiplex", "genome_alignment")
  ppl@durations <- d
  ppl@last_error <- list(step = "gene_quantification", error = "boom")
  # an existing output file exercises the display_outputs "exists" + size branch
  bam <- file.path(ppl@outdir, "genome.bam")
  writeLines("x", bam)
  ppl@genome_bam <- bam

  out <- cli::cli_fmt(show(ppl))
  expect_true(any(grepl("completed", out)))
  expect_true(any(grepl("failed", out)))
  expect_true(any(grepl("pending", out)))
})

test_that("show() handles a MultiSample object with vector output slots", {
  ppl <- example_pipeline("MultiSampleSCPipeline")
  # mixed existing / missing per-sample outputs exercises the multi-path branch
  existing <- file.path(ppl@outdir, "s1.bam")
  writeLines("x", existing)
  ppl@genome_bam <- c(existing, file.path(ppl@outdir, "missing.bam"))
  expect_no_error(show(ppl))
})
