# Additional config_functions.R branch coverage: load_config type dispatch,
# deprecated-pattern warning, unknown-parameter error, and the remaining
# check_arguments validation branches.

test_that("load_config supports SIRV defaults and rejects unknown types", {
  outdir <- withr::local_tempdir()
  cfg <- file.path(outdir, "min.json")
  write(jsonlite::toJSON(list(pipeline_parameters = list(threads = 4)), pretty = TRUE), cfg)

  sirv <- load_config(cfg, type = "SIRV")
  expect_true(sirv$alignment_parameters$no_flank)
  expect_error(load_config(cfg, type = "invalid"), "Unrecognised config type")
})

test_that("merge_configs_recursive warns on the deprecated pattern field and preserves it", {
  default <- jsonlite::fromJSON(
    system.file("extdata", "config_sclr_nanopore_3end.json", package = "FLAMES")
  )
  user <- list(barcode_parameters = list(pattern = list(BC = "NNNN")))
  expect_warning(
    merged <- FLAMES:::merge_configs_recursive(default, user),
    "deprecated"
  )
  expect_equal(merged$barcode_parameters$pattern, list(BC = "NNNN"))
})

test_that("create_config errors on an unknown parameter", {
  expect_error(
    create_config(withr::local_tempdir(), totally_bogus_param = 1),
    "Unknown parameter"
  )
})

test_that("check_arguments creates outdir / defaults config and validates inputs", {
  fa <- system.file("extdata", "rps24.fa.gz", package = "FLAMES")
  gtf <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")

  # non-existent outdir is created; NULL config_file -> a default is generated
  outdir <- file.path(withr::local_tempdir(), "new_sub")
  res <- FLAMES:::check_arguments(
    annotation = gtf, fastq = NULL, genome_bam = NULL,
    outdir = outdir, genome_fa = fa, config_file = NULL
  )
  expect_true(dir.exists(outdir))
  expect_type(res, "list")

  cfg <- create_config(withr::local_tempdir())
  expect_error(
    FLAMES:::check_arguments(gtf, "/no/such.fastq", NULL, withr::local_tempdir(), fa, cfg),
    "exists"
  )
  expect_error(
    FLAMES:::check_arguments(gtf, NULL, "/no/such.bam", withr::local_tempdir(), fa, cfg),
    "genome_bam"
  )
  cfg_bad <- create_config(withr::local_tempdir(), isoform_parameters.downsample_ratio = 2)
  expect_error(
    FLAMES:::check_arguments(gtf, NULL, NULL, withr::local_tempdir(), fa, cfg_bad),
    "downsample_ratio"
  )
})
