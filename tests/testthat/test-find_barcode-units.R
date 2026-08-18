# tests for find_barcode.R helpers (not the compiled stuff)

# ---- barcode_segment -------------------------------------------------------

test_that("barcode_segment validates arguments and builds a MATCHED_SPLIT segment", {
  expect_error(barcode_segment("BADTYPE", "NNNN", "x"), "type must be one of")
  expect_error(barcode_segment("FIXED", pattern = ""), "non-empty character")
  expect_error(barcode_segment("MATCHED", "NNNN", name = ""), "name must be")
  expect_error(barcode_segment("MATCHED", "NNNN", "CB", bc_list = "CB", buffer_size = -1), "non-negative")
  expect_error(barcode_segment("MATCHED", "NNNN", "CB", bc_list = "CB", max_edit_distance = -1), "non-negative")
  expect_error(barcode_segment("MATCHED", "NNNN", "CB"), "bc_list must be provided")
  expect_error(barcode_segment("MATCHED_SPLIT", "NNNN", "CB"), "group must be provided")
  expect_error(barcode_segment("MATCHED_SPLIT", "NNNN", "CB", group = "g", bc_list = "other"),
               "same as group")

  seg <- barcode_segment("MATCHED_SPLIT", "NNNN", "CB", group = "g")
  expect_s4_class(seg, "FlexiplexSegment")
  expect_equal(seg@bc_list_name, "g") # bc_list defaults to the group name
  expect_equal(barcode_segment("FIXED", "ACGT")@name, "FIXED_SEGMENT")
})

# ---- barcode_group ---------------------------------------------------------

test_that("barcode_group validates arguments and builds a group", {
  expect_error(barcode_group("", "CB"), "name must be")
  expect_error(barcode_group("g", ""), "bc_list_name must be")
  expect_error(barcode_group("g", "CB", max_edit_distance = "x"), "numeric scalar")
  expect_error(barcode_group("g", "CB", max_edit_distance = -1), "non-negative")

  grp <- barcode_group("g", "CB")
  expect_s4_class(grp, "FlexiplexGroup")
  expect_equal(grp@name, "g")
})

# ---- .legacy_pattern_to_flexiplex_segments ---------------------------------

test_that(".legacy_pattern_to_flexiplex_segments validates the pattern", {
  f <- FLAMES:::.legacy_pattern_to_flexiplex_segments
  expect_error(f(42, "file", 2, 2), "named character")
  expect_error(f(c("ACGT"), "file", 2, 2), "named character")
  expect_error(f(c(primer = "ACGT"), "file", 2, 2), "'BC'")
  expect_error(f(c(BC = "NNNN"), "", 2, 2), "single non-empty")
})

# ---- .resolve_bc_lists -----------------------------------------------------

test_that(".resolve_bc_lists resolves named keys, file paths, groups, and errors", {
  seg <- barcode_segment("MATCHED", "NNNN", "CB", bc_list = "CB")

  expect_error(FLAMES:::.resolve_bc_lists(list(seg), 42), "character vector")
  expect_error(FLAMES:::.resolve_bc_lists(list(seg), c("a", "b")), "named character vector")
  expect_error(FLAMES:::.resolve_bc_lists(list("nope"), c(CB = "/x")),
               "FlexiplexSegment or FlexiplexGroup")

  # named lookup
  r <- FLAMES:::.resolve_bc_lists(list(seg), c(CB = "/path/to/bc"))
  expect_equal(r[[1]]@bc_list_name, "/path/to/bc")

  # bc_list_name that is itself an existing file -> used directly (with a message)
  real <- system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES")
  seg_file <- barcode_segment("MATCHED", "NNNN", "CB", bc_list = real)
  expect_message(
    r_file <- FLAMES:::.resolve_bc_lists(list(seg_file), c(OTHER = "/x")),
    "using it as a file path"
  )
  expect_equal(r_file[[1]]@bc_list_name, real)

  # unresolvable key
  seg_missing <- barcode_segment("MATCHED", "NNNN", "CB", bc_list = "MISSING")
  expect_error(FLAMES:::.resolve_bc_lists(list(seg_missing), c(OTHER = "/x")),
               "does not match any name")

  # FlexiplexGroup branch
  grp <- barcode_group("g", "CB")
  rg <- FLAMES:::.resolve_bc_lists(list(grp), c(CB = "/path/g"))
  expect_equal(rg[[1]]@bc_list_name, "/path/g")
})

# ---- find_barcode strand guard ---------------------------------------------

test_that("find_barcode rejects an invalid strand before doing any work", {
  expect_error(
    find_barcode(
      fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
      barcodes_files = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
      strand = "?"
    ),
    "strand must be"
  )
})

# ---- plot_demultiplex_raw --------------------------------------------------

# a small in-memory find_barcode result (one "sample")
make_demux_result <- function(with_cutadapt = FALSE) {
  res <- list(
    read_counts = c("total reads" = 100, "demultiplexed reads" = 80, "single match reads" = 70),
    reads_tb = data.frame(
      CB = rep(c("AAAA", "CCCC", "GGGG"), c(5, 3, 2)),
      UB = paste0("U", 1:10),
      FlankEditDist = c(0, 0, 1, 1, 2, 0, 1, 2, 0, 1)
    )
  )
  if (with_cutadapt) {
    res$cutadapt <- list(read_counts = list(input = 100, output = 90, read1_with_adapter = 88))
  }
  res
}

test_that("plot_demultiplex_raw returns ggplots for a single in-memory sample", {
  p <- plot_demultiplex_raw(make_demux_result())
  expect_type(p, "list")
  expect_s3_class(p$reads_count_plot, "ggplot")
  expect_s3_class(p$knee_plot, "ggplot")
  expect_s3_class(p$flank_editdistance_plot, "ggplot")
  expect_null(p$cutadapt_plot)
})

test_that("plot_demultiplex_raw handles multiple samples and a cutadapt report", {
  res <- list(S1 = make_demux_result(TRUE), S2 = make_demux_result(TRUE))
  p <- plot_demultiplex_raw(res)
  expect_s3_class(p$knee_plot, "ggplot")
  expect_s3_class(p$cutadapt_plot, "ggplot")
})

test_that("plot_demultiplex_raw reads reads_tb from stats_out when absent", {
  res <- list(
    read_counts = c("total reads" = 100, "demultiplexed reads" = 80, "single match reads" = 70),
    stats_out = test_path("bc_stat")
  )
  p <- suppressWarnings(plot_demultiplex_raw(res))
  expect_s3_class(p$knee_plot, "ggplot")
})

test_that("plot_demultiplex method errors before the demultiplex step has run", {
  ppl <- example_pipeline("SingleCellPipeline")
  expect_error(plot_demultiplex(ppl), "have you run the demultiplex step")
})
