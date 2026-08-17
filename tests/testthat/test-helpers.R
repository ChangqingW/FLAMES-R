# -- minimap2_align.R

test_that("check_status_code passes on 0 and errors on a non-zero status", {
  expect_no_error(FLAMES:::check_status_code(0, "some command"))
  expect_error(
    suppressWarnings(FLAMES:::check_status_code(1, "some command", "MyProc")),
    "status code 1"
  )
})

test_that("gff2bed writes a BED file from a GTF", {
  out <- file.path(withr::local_tempdir(), "out.bed")
  gff <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
  FLAMES:::gff2bed(gff, out)
  expect_true(file.exists(out))
  expect_gt(length(readLines(out)), 0)
})

test_that("plot_flagstat returns a ggplot", {
  fs <- data.frame(total = 10, mapped = 8, primary = 8, secondary = 0)
  expect_s3_class(FLAMES:::plot_flagstat(fs), "ggplot")
})

test_that("get_flagstat parses flag summary via Rsamtools when samtools is absent", {
  bam <- local_bam()
  df <- FLAMES:::get_flagstat(bam, NA)
  expect_true(is.data.frame(df))
  expect_true(all(c("total", "mapped", "primary", "secondary") %in% names(df)))
  expect_equal(df[["total"]], 4)   # local_bam() writes 4 reads
  expect_equal(df[["mapped"]], 4)  # all mapped (flag 0)
})

# -- quantification.R

test_that("run_oarfish returns the output stem on success and stops on failure", {
  bam <- local_bam()
  dir <- withr::local_tempdir()
  # `echo` exits 0: success; `false` exits 1: stop
  out <- FLAMES:::run_oarfish(bam, dir, oarfish_bin = "echo",
                              single_cell = TRUE, additional_args = character(0))
  expect_equal(out, file.path(dir, "oarfish"))
  expect_error(
    FLAMES:::run_oarfish(bam, dir, oarfish_bin = "false", additional_args = character(0)),
    "error running oarfish"
  )
})

test_that("quantify_gene validates BAM inputs before any python call", {
  gtf <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
  dir <- withr::local_tempdir()
  expect_error(
    FLAMES:::quantify_gene(gtf, dir, in_bam = file.path(dir, "nope.bam")),
    "missing"
  )
  bam1 <- local_bam(); bam2 <- local_bam()
  expect_error(
    FLAMES:::quantify_gene(gtf, dir, pipeline = "sc_single_sample", in_bam = c(bam1, bam2)),
    "Incorrect number"
  )
})

test_that("quantify_transcript_flames stops when no realignment bam is present", {
  gtf <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
  dir <- withr::local_tempdir() # no *realign2transcript.bam
  expect_error(
    FLAMES:::quantify_transcript_flames(gtf, dir, config = list(), pipeline = "sc_single_sample"),
    "Incorrect number of realignment"
  )
})

# -- find_isoform.R

test_that("find_isoform_flames stops when the bam index is missing", {
  expect_error(
    FLAMES:::find_isoform_flames("ann.gtf", "genome.fa", "/nonexistent/x.bam", "outdir", list()),
    "Cannot find corresponding"
  )
})
