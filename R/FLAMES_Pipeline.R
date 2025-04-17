setClass(
  "FLAMES.SingleCellPipeline",
  contains = "FLAMES.Pipeline",
  slots = list(
    # inputs
    barcodes_file = "character",        # Path to the barcodes file
    expect_cell_number = "numeric",     # Expected number of cells
    # outputs
    demultiplexed_fastq = "character",  # path to demultiplexed FASTQ files
    dedputed_fastq = "character",       # path to deduplicated FASTQ files
    metadata = "list"                   # Metadata for the pipeline run (metrics, etc.)
  ),
  prototype = list(
    barcodes_file = NA_character_,
    expect_cell_number = NA_real_,
    demultiplexed_fastq = NA_character_,
    dedputed_fastq = NA_character_,
    metadata = list()
  )
)

SingleCellPipeline <- function(
  config_file, outdir, fastq, annotation, genome_fa, barcodes_file, expect_cell_number
) {
  pipeline <- new("FLAMES.SingleCellPipeline")
  config <- check_arguments(annotation, fastq, genome_bam = NULL, outdir, genome_fa, config_file)$config

  if (!dir.exists(outdir)) {
    dir.create(outdir)
    message(sprintf("Output directory (%s) did not exist, created.", outdir))
  }

  steps <- c(
    "barcode_demultiplex", "gene_quantification", "genome_alignment",
    "isoform_identification", "read_realignment", "transcript_quantification"
  )
  steps <- config$pipeline_parameters[paste0("do_", steps)] |>
    unlist() |>
    setNames(steps)
  message("Configured steps: \n", paste0("\t", names(steps), ": ", steps, collapse = "\n"))

  # assign slots
  ## inputs
  pipeline@config <- config
  pipeline@outdir <- outdir
  pipeline@fastq <- fastq
  pipeline@annotation <- annotation
  pipeline@genome_fa <- genome_fa
  if (!missing(barcodes_file)) {
    pipeline@barcodes_file <- barcodes_file
  } else if (!missing(expect_cell_number)) {
    pipeline@expect_cell_number <- expect_cell_number
  } else {
    stop("Either barcodes_file or expect_cell_number must be provided.")
  }
  ## outputs
  # metadata
  pipeline@genome_bam <- file.path(outdir, "align2genome.bam")
  pipeline@transcript_bam <- file.path(outdir, "align2transcript.bam")
  ## pipeline state
  pipeline@steps <- steps
  pipeline@completed_steps <- setNames(
    rep(FALSE, length(steps)), names(steps)
  )
  pipeline@demultiplexed_fastq <- file.path(outdir, "matched_reads.fastq")
  pipeline@dedputed_fastq <- file.path(outdir, "matched_reads_dedup.fastq")

  # TODO: add resume option
  # validate if e.g. genome_bam exists, skip genome alignment step
  return(pipeline)
}

# Multi-sample pipeline is the same as SingleCellPipeline
# but input slots will be vectors of the same length
setClass(
  "FLAMES.MultiSampleSCPipeline",
  contains = "FLAMES.SingleCellPipeline"
)

MultiSampleSCPipeline <- function(
  config_file, outdir, fastq, annotation, genome_fa, barcodes_file, expect_cell_number
) {
  pipeline <- new("FLAMES.MultiSampleSCPipeline")
  config <- check_arguments(annotation, fastq, genome_bam = NULL, outdir, genome_fa, config_file)$config

  if (!dir.exists(outdir)) {
    dir.create(outdir)
    message(sprintf("Output directory (%s) did not exist, created.", outdir))
  }

  # If fastq is a directory, replace it with fastq files in the directory
  if (length(fastq) == 1) {
    if (utils::file_test("-f", fastq)) {
      stop("Only one fastq file provided, but multiple samples are expected.")
    }
    fastq <- list.files(fastq, pattern = "\\.(fastq|fq)(\\.gz)?$", full.names = TRUE)
  }
  # Use basename if fastq is unnamed
  if (is.null(names(fastq))) {
    names(fastq) <- gsub("\\.(fastq|fq)(\\.gz)?$", "", basename(fastq))
  }

  # Check if barcodes_file matches fastq
  if (!missing(barcodes_file)) {
    if (length(barcodes_file) == 1) {
      message("Only one barcodes file provided, assuming it applies to all samples.")
      barcodes_file <- rep(barcodes_file, length(fastq)) |>
        setNames(names(fastq))
    } else if (length(barcodes_file) != length(fastq)) {
      stop("Number of barcodes files must match the number of fastq files.")
    } else if (is.null(names(barcodes_file))) {
      message("Assuming barcodes files are provided in the same order as the corresponding fastq files.")
      names(barcodes_file) <- names(fastq)
    } else if (any(!names(barcodes_file) %in% names(fastq))) {
      stop(
        sprintf(
          "Barcodes files names (%s) do not match fastq files names (%s).",
          paste(names(barcodes_file), collapse = ", "),
          paste(names(fastq), collapse = ", ")
        )
      )
    }
  }
  # same but for expect_cell_number
  if (!missing(expect_cell_number)) {
    if (length(expect_cell_number) == 1) {
      expect_cell_number <- rep(expect_cell_number, length(fastq)) |>
        setNames(names(fastq))
    } else if (length(expect_cell_number) != length(fastq)) {
      stop("expect_cell_number must be a single value (applied to all sample) or a vector of the same length as fastq inputs.")
    } else if (is.null(names(expect_cell_number))) {
      names(expect_cell_number) <- names(fastq)
    } else if (any(!names(expect_cell_number) %in% names(fastq))) {
      stop(
        sprintf(
          "names(expect_cell_number) (%s) do not match fastq files names (%s).",
          paste(names(expect_cell_number), collapse = ", "),
          paste(names(fastq), collapse = ", ")
        )
      )
    }
  }

  steps <- c(
    "barcode_demultiplex", "gene_quantification", "genome_alignment",
    "isoform_identification", "read_realignment", "transcript_quantification"
  )
  steps <- config$pipeline_parameters[paste0("do_", steps)] |>
    unlist() |>
    setNames(steps)
  message("Configured steps: \n", paste0("\t", names(steps), ": ", steps, collapse = "\n"))

  # assign slots
  ## inputs
  pipeline@config <- config
  pipeline@outdir <- outdir
  pipeline@fastq <- fastq
  pipeline@annotation <- annotation
  pipeline@genome_fa <- genome_fa
  if (!missing(barcodes_file)) {
    pipeline@barcodes_file <- barcodes_file
  } else if (!missing(expect_cell_number)) {
    pipeline@expect_cell_number <- expect_cell_number
  } else {
    stop("Either barcodes_file or expect_cell_number must be provided.")
  }
  ## outputs
  # metadata
  pipeline@genome_bam <- file.path(outdir, paste0(names(fastq), "_", "align2genome.bam"))
  pipeline@transcript_bam <- file.path(outdir, paste0(names(fastq), "_", "align2transcript.bam"))
  ## pipeline state
  pipeline@steps <- steps
  pipeline@completed_steps <- setNames(
    rep(FALSE, length(steps)), names(steps)
  )
  pipeline@demultiplexed_fastq <- file.path(outdir, paste0(names(fastq), "_matched_reads.fastq"))
  pipeline@dedputed_fastq <- file.path(outdir, paste0(names(fastq), "_matched_reads_dedup.fastq"))

  # TODO: add resume option
  return(pipeline)
}

