# getters/setters, validation, controllers, run/resume/plot etc.

# -- getters/setters

test_that("steps() getter/setter validate names", {
  ppl <- example_pipeline("BulkPipeline")
  expect_type(steps(ppl), "logical")
  expect_true(all(c("genome_alignment", "transcript_quantification") %in% names(steps(ppl))))

  steps(ppl)["genome_alignment"] <- FALSE
  expect_false(steps(ppl)[["genome_alignment"]])

  expect_error({ steps(ppl) <- c(TRUE, FALSE) }, "named logical vector")
  expect_error({ steps(ppl) <- c(not_a_step = TRUE) }, "Invalid step names")
})

test_that("config() getter/setter accepts a list, a file path, and rejects other input", {
  ppl <- example_pipeline("BulkPipeline")
  expect_type(config(ppl), "list")
  expect_true("pipeline_parameters" %in% names(config(ppl)))

  cfg_path <- create_config(withr::local_tempdir(), threads = 7)
  config(ppl) <- cfg_path
  expect_equal(config(ppl)$pipeline_parameters$threads, 7)

  config(ppl) <- list(pipeline_parameters = list(threads = 3))
  expect_equal(config(ppl)$pipeline_parameters$threads, 3)

  expect_error({ config(ppl) <- 42 }, "must be a list or a path")
})

test_that("experiment() returns NULL when unset/missing and the object when present", {
  ppl <- example_pipeline("BulkPipeline")
  expect_null(experiment(ppl))

  ppl@experiment <- file.path(ppl@outdir, "nope.rds")
  expect_warning(res <- experiment(ppl), "not found")
  expect_null(res)

  f <- file.path(ppl@outdir, "exp.rds")
  saveRDS(1:5, f)
  ppl@experiment <- f
  expect_equal(experiment(ppl), 1:5)
})

# -- controllers

test_that("normalize_controllers wraps, warns on unknown names, and rejects bad input", {
  skip_if_not_installed("crew")
  ctrl <- crew::crew_controller_local()

  wrapped <- FLAMES:::normalize_controllers(ctrl, c("genome_alignment", "read_realignment"))
  expect_named(wrapped, "default")

  expect_warning(
    FLAMES:::normalize_controllers(list(bogus_step = ctrl), c("genome_alignment")),
    "do not match"
  )
  expect_error(FLAMES:::normalize_controllers("nope", c("genome_alignment")), "crew_class_controller")
})

test_that("controllers() getter/setter round-trips and validates", {
  skip_if_not_installed("crew")
  ppl <- example_pipeline("BulkPipeline")
  expect_type(controllers(ppl), "list")

  controllers(ppl) <- crew::crew_controller_local()
  expect_named(controllers(ppl), "default")

  expect_error({ controllers(ppl) <- "nope" }, "crew_class_controller")
})

# -- prerun_check

test_that("prerun_check handles fresh / completed / partial states", {
  ppl <- example_pipeline("BulkPipeline")
  expect_true(FLAMES:::prerun_check(ppl))

  done <- ppl
  done@completed_steps[] <- TRUE
  expect_message(expect_false(FLAMES:::prerun_check(done, overwrite = FALSE)), "already been completed")
  expect_warning(expect_true(FLAMES:::prerun_check(done, overwrite = TRUE)), "Re-running")

  partial <- ppl
  partial@completed_steps[1] <- TRUE
  expect_error(FLAMES:::prerun_check(partial, overwrite = FALSE), "partially completed")
  expect_warning(expect_true(FLAMES:::prerun_check(partial, overwrite = TRUE)), "Re-running")
})

# -- run_FLAMES / resume_FLAMES early-return branches

test_that("run_FLAMES returns early when the pipeline is already complete", {
  ppl <- example_pipeline("BulkPipeline")
  ppl@completed_steps[] <- TRUE
  expect_message(res <- run_FLAMES(ppl, overwrite = FALSE), "already been completed")
  expect_s4_class(res, "FLAMES.Pipeline")
})

test_that("resume_FLAMES messages and returns when nothing is left to do", {
  ppl <- example_pipeline("BulkPipeline")
  ppl@completed_steps[] <- TRUE
  expect_message(res <- resume_FLAMES(ppl), "already been completed")
  expect_s4_class(res, "FLAMES.Pipeline")
})

# -- run_step dispatch: unknown step / base-class dummy stops

test_that("run_step errors on an unknown step and on unimplemented bulk steps", {
  ppl <- example_pipeline("BulkPipeline")
  expect_error(run_step(ppl, "not_a_real_step"), "Unknown step")
  expect_error(run_step(ppl, "barcode_demultiplex"), "single cell")
  expect_error(run_step(ppl, "gene_quantification"), "not implemented for bulk")
})

# -- plot_durations

test_that("plot_durations errors without durations and returns a ggplot with them", {
  ppl <- example_pipeline("BulkPipeline")
  expect_error(plot_durations(ppl), "No step durations")

  d <- as.difftime(c(5, 130, 12), units = "secs")
  names(d) <- c("genome_alignment", "isoform_identification", "read_realignment")
  ppl@durations <- d
  expect_s3_class(plot_durations(ppl), "ggplot")
})

# -- check_arguments validation

test_that("check_arguments errors on a missing annotation file", {
  outdir <- withr::local_tempdir()
  genome_fa <- system.file("extdata", "rps24.fa.gz", package = "FLAMES")
  cfg <- create_config(outdir)
  expect_error(
    FLAMES:::check_arguments(
      annotation = file.path(outdir, "nonexistent.gtf"), fastq = NULL,
      genome_bam = NULL, outdir = outdir, genome_fa = genome_fa, config_file = cfg
    ),
    "exists"
  )
})

test_that("check_arguments enforces GTF for bambu and warns on oarfish without gene quant", {
  outdir <- withr::local_tempdir()
  gtf <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
  fa <- system.file("extdata", "rps24.fa.gz", package = "FLAMES")

  # bambu requires a GTF: give an (existing) non-GTF annotation
  cfg_bambu <- create_config(outdir, pipeline_parameters.bambu_isoform_identification = TRUE)
  expect_error(
    FLAMES:::check_arguments(
      annotation = fa, fastq = NULL, genome_bam = NULL,
      outdir = outdir, genome_fa = fa, config_file = cfg_bambu
    ),
    "Bambu requires GTF"
  )

  cfg_oarfish <- create_config(
    outdir,
    pipeline_parameters.bambu_isoform_identification = FALSE,
    pipeline_parameters.do_transcript_quantification = TRUE,
    pipeline_parameters.oarfish_quantification = TRUE,
    pipeline_parameters.do_gene_quantification = FALSE
  )
  expect_warning(
    FLAMES:::check_arguments(
      annotation = gtf, fastq = NULL, genome_bam = NULL,
      outdir = outdir, genome_fa = fa, config_file = cfg_oarfish
    ),
    "oarfish"
  )
})

multi_inputs <- function(dir, n = 3) {
  fq_dir <- file.path(dir, "fastq")
  dir.create(fq_dir, showWarnings = FALSE)
  src <- system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES")
  fqs <- file.path(fq_dir, paste0("sample", seq_len(n), ".fastq.gz"))
  for (f in fqs) file.copy(src, f, overwrite = TRUE)
  genome_fa <- file.path(dir, "rps24.fa")
  R.utils::gunzip(system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
                  destname = genome_fa, remove = FALSE, overwrite = TRUE)
  list(
    fastq = fqs, genome_fa = genome_fa,
    annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
    config = create_config(dir), outdir = dir
  )
}

build_multi <- function(io, ...) {
  MultiSampleSCPipeline(
    config_file = io$config, outdir = io$outdir, fastq = io$fastq,
    annotation = io$annotation, genome_fa = io$genome_fa,
    minimap2 = "minimap2", samtools = "samtools", ...
  )
}

test_that("MultiSampleSCPipeline errors on barcodes_file / expect_cell_number length mismatch", {
  io <- multi_inputs(withr::local_tempdir(), 3)
  expect_error(build_multi(io, barcodes_file = c("a", "b")), "does not match number of fastq")
  expect_error(build_multi(io, expect_cell_number = c(100, 200)), "does not match number of fastq")
})

test_that("MultiSampleSCPipeline rejects a single fastq file", {
  io <- multi_inputs(withr::local_tempdir(), 1)
  expect_error(build_multi(io), "single-sample")
})

test_that("MultiSampleSCPipeline errors when the fastq folder has too few files", {
  dir <- withr::local_tempdir()
  fq_dir <- file.path(dir, "fq")
  dir.create(fq_dir)
  file.copy(system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
            file.path(fq_dir, "only.fastq.gz"))
  genome_fa <- file.path(dir, "rps24.fa")
  R.utils::gunzip(system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
                  destname = genome_fa, remove = FALSE, overwrite = TRUE)
  expect_error(
    MultiSampleSCPipeline(
      config_file = create_config(dir), outdir = dir, fastq = fq_dir,
      annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
      genome_fa = genome_fa, minimap2 = "minimap2", samtools = "samtools"
    ),
    "file\\(s\\) found"
  )
})

test_that("MultiSampleSCPipeline requires expect_cell_number when demultiplexer is BLAZE", {
  dir <- withr::local_tempdir()
  io <- multi_inputs(dir, 3)
  cfg <- create_config(dir, pipeline_parameters.demultiplexer = "BLAZE")
  expect_error(
    MultiSampleSCPipeline(
      config_file = cfg, outdir = io$outdir, fastq = io$fastq,
      annotation = io$annotation, genome_fa = io$genome_fa,
      minimap2 = "minimap2", samtools = "samtools"
    ),
    "expect_cell_number"
  )
})

test_that("example_pipeline rejects an unknown pipeline type", {
  expect_error(example_pipeline("not_a_pipeline"), "Invalid pipeline type")
})

# -- step-method demultiplex parameter validation (errors before any tool runs)

test_that("barcode_demultiplex errors on an unknown demultiplexer", {
  for (type in c("SingleCellPipeline", "MultiSampleSCPipeline")) {
    ppl <- example_pipeline(type)
    ppl@config$pipeline_parameters$demultiplexer <- "nonsense"
    ppl@controllers <- list()
    expect_error(run_step(ppl, "barcode_demultiplex"), "[Uu]nknown demultiplexer")
  }
})

test_that("barcode_demultiplex requires expect_cell_number for BLAZE", {
  for (type in c("SingleCellPipeline", "MultiSampleSCPipeline")) {
    ppl <- example_pipeline(type)
    ppl@config$pipeline_parameters$demultiplexer <- "BLAZE"
    ppl@expect_cell_number <- NA_real_
    ppl@controllers <- list()
    expect_error(run_step(ppl, "barcode_demultiplex"), "expect_cell_number")
  }
})
