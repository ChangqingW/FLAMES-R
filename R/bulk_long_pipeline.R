# Define the base S4 class for FLAMES pipelines
setClass(
  "FLAMES.Pipeline",
  slots = list(
    # inputs
    config = "list",            # Configuration parameters
    outdir = "character",       # Output directory
    fastq = "character",        # Path to the FASTQ files
    annotation = "character",   # Path to the annotation file
    genome_fa = "character",    # Path to the genome FASTA file

    # outputs
    metadata = "list",          # Metadata for the pipeline
    genome_bam = "character",     # Path to the genome BAM file
    transcriptome_bam = "character", # Path to the transcript BAM file
    novel_isoform_annotation = "character", # Path to the novel isoform GFF / GTF file
    transcriptome_assembly = "character", # Path to the transcriptome assembly file
    experiment = "SummarizedExperiment", # SummarizedExperiment object for quantification results

    # binaries
    minimap2 = "character",     # Path to the minimap2 binary
    k8 = "character",           # Path to the k8 binary
    samtools = "character",     # Path to the samtools binary

    # pipeline state
    steps = "character",        # Steps to perform
    completed_steps = "character" # Completed steps
  ),
  prototype = list(
    config = list(),
    outdir = NA_character_,
    fastq = NA_character_,
    annotation = NA_character_,
    genome_fa = NA_character_,

    metadata = list(),
    genome_bam = NA_character_,
    transcriptome_bam = NA_character_,
    novel_isoform_annotation = NA_character_,
    transcriptome_assembly = NA_character_,
    experiment = NA,

    minimap2 = NA_character_,
    k8 = NA_character_,
    samtools = NA_character_,

    steps = character(),
    completed_steps = character()
  )
)

# constructor for the bulk pipeline
BulkPipeline <- function(config_file, outdir, fastq, annotation, genome_fa, minimap2, samtools, k8) {

  pipeline <- new("FLAMES.Pipeline")
  config <- check_arguments(annotation, fastq, genome_bam = NULL, outdir, genome_fa, config_file)$config

  if (!dir.exists(outdir)) {
    dir.create(outdir)
    message(sprintf("Output directory (%s) did not exist, created.", outdir))
  }

  if (utils::file_test("-d", fastq)) {
    fastq <- list.files(fastq, pattern = "\\.(fastq|fq)(\\.gz)?$", full.names = TRUE)
  } else if (!utils::file_test("-f", fastq)) {
    stop("fastq must be a valid path to a folder or a FASTQ file")
  }
  names(fastq) <- gsub("\\.(fastq|fq)(\\.gz)?$", "", basename(fastq))

  steps <- c(
    "genome_alignment", "isoform_identification",
    "read_realignment", "transcript_quantification"
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

  ## outputs
  # metadata
  pipeline@genome_bam <- file.path(outdir, paste0(names(fastq), "_", "align2genome.bam"))
  pipeline@transcriptome_bam <- file.path(outdir, paste0(names(fastq), "_", "align2transcript.bam"))
  pipeline@transcriptome_assembly <- file.path(outdir, "transcriptome_assembly.fa")

  ## binaries
  if (missing(minimap2)) {
    minimap2 <- find_bin("minimap2")
    if (is.na(minimap2)) {
      stop("minimap2 not found, please make sure it is installed and provide its path as the minimap2 argument")
    }
  }
  pipeline@minimap2 <- minimap2
  if (missing(k8)) {
    k8 <- find_bin("k8")
    if (is.na(k8)) {
      stop("k8 not found, please make sure it is installed and provide its path as the k8 argument")
    }
  }
  pipeline@k8 <- k8
  if (missing(samtools)) {
    samtools <- find_bin("samtools")
  }
  if (is.na(samtools)) {
    message("samtools not found, will use Rsamtools package instead")
  }
  pipeline@samtools <- samtools

  ## pipeline state
  pipeline@steps <- steps
  pipeline@completed_steps <- setNames(
    rep(FALSE, length(steps)), names(steps)
  )
  # TODO: add resume option
  # validate if e.g. genome_bam exists, skip genome alignment step

  return(pipeline)
}

setGeneric("prerun_check", function(pipeline) {
  standardGeneric("prerun_check")
})
setMethod("prerun_check", "FLAMES.Pipeline", function(pipeline, overwrite = FALSE) {
  # return TRUE for run_FLAMES to proceed
  # return FALSE for run_FLAMES to stop gracefully, when all steps are completed
  # stop() when pipeline is partially completed and user does not want to overwrite
  if (all(pipeline@completed_steps)) {
    message("All steps have already been completed.")
    if (overwrite) {
      warning("Re-running pipeline to overwrite existing results.")
      return(TRUE)
    } else {
      message("Pipeline is already completed. Set overwrite = TRUE to re-run.")
      return(FALSE)
    }
  } else if (any(pipeline@completed_steps)) {
    message("Some steps have already been completed.")
    if (overwrite) {
      warning("Re-running pipeline to overwrite existing results.")
      return(TRUE)
    } else {
      # TODO: implement resuming and prompt user to resume
      stop("Pipeline is partially completed. Please set overwrite = TRUE to proceed.")
    }
  } else {
    # nothing done yet, proceed
    return(TRUE)
  }
})

setGeneric("run_step", function(pipeline, step) {
  standardGeneric("run_step")
})
setMethod("run_step", "FLAMES.Pipeline", function(pipeline, step) {
  pipeline <- switch(step,
    barcode_demultiplexing = barcode_demultiplexing(pipeline),
    genome_alignment = genome_alignment(pipeline),
    isoform_identification = isoform_identification(pipeline),
    read_realignment = read_realignment(pipeline),
    transcript_quantification = transcript_quantification(pipeline),
    stop(sprintf("Unknown step: %s", step))
  )
  pipeline@completed_steps[step] <- TRUE
  return(pipeline)
})

# individual steps as methods
setGeneric("barcode_demultiplexing", function(pipeline) {
  standardGeneric("barcode_demultiplexing")
})
setMethod("barcode_demultiplexing", "FLAMES.Pipeline", function(pipeline) {
  stop("Barcode demultiplexing is only implemented for single cell pipelines (SingleCellPipeline() and MultiSampleSCPipeline())")
})
setGeneric("genome_alignment", function(pipeline) {
  standardGeneric("genome_alignment")
})
setMethod("genome_alignment", "FLAMES.Pipeline", function(pipeline) {

  minimap2_args <- c(
    "-ax", "splice",  "-k14", "--secondary=no", "-y",
    "-t", pipeline@config$pipeline_parameters$threads,
    "--seed", pipeline@config$pipeline_parameters$seed
  )
  if (pipeline@config$alignment_parameters$no_flank) {
    minimap2_args <- base::append(minimap2_args, "--splice-flank=no")
  }

  # k8 paftools.js gff2bend gff > bed12
  paftoolsjs <- system.file("paftools.js", package = "FLAMES")
  if (pipeline@config$alignment_parameters$use_junctions) {
    bed_file <- tempfile(tmpdir = pipeline@outdir, fileext = ".bed")
    paftoolsjs_status <- base::system2(
      command = pipeline@k8,
      args = c(paftoolsjs, "gff2bed", pipeline@annotation, ">", bed_file)
    )
    if (!is.null(base::attr(paftoolsjs_status, "status")) && base::attr(paftoolsjs_status, "status") != 0) {
      stop(sprintf(
        "Error running %s\nAre you using NCBI GFF3? It is not well supported by minimap2's paftools.js, see https://github.com/lh3/minimap2/issues/422",
        paste(c(pipeline@k8, paftoolsjs, "gff2bed", pipeline@annotation, ">", bed_file), collapse = " ")
      ))
    }
    minimap2_args <- base::append(minimap2_args, c("--junc-bed", bed_file, "--junc-bonus", "1"))
  }

  res <- lapply(
    seq_along(pipeline@fastq),
    function(i) {
      if (!is.null(names(pipeline@fastq))) {
        sample <- names(pipeline@fastq)[i]
      } else {
        sample <- pipeline@fastq[i]
      }
      message(sprintf("Aligning sample %s -> %s", sample, pipeline@genome_bam[i]))
      minimap2_align(
        fq_in = pipeline@fastq[i],
        fa_file = pipeline@genome_fa,
        config = pipeline@config,
        outfile = pipeline@genome_bam[i],
        minimap2_args = minimap2_args,
        minimap2 = pipeline@minimap2,
        samtools = pipeline@samtools,
        threads = pipeline@config$pipeline_parameters$threads
      )
    }
  )
  if (!is.null(names(pipeline@fastq))) {
    names(res) <- names(pipeline@fastq)
  }
  unlink(bed_file)
  pipeline@metadata$genome_alignment <- res
  return(pipeline)
})

setGeneric("isoform_identification", function(pipeline) {
  standardGeneric("isoform_identification")
})
setMethod("isoform_identification", "FLAMES.Pipeline", function(pipeline) {
  if (pipeline@config$pipeline_parameters$bambu_isoform_identification) {
    novel_isoform_annotation <- find_isoform_bambu(
      annotation = pipeline@annotation,
      genome_fa = pipeline@genome_fa,
      genome_bam = pipeline@genome_bam,
      outdir = pipeline@outdir,
      config = pipeline@config
    )
  } else {
    novel_isoform_annotation <- find_isoform_flames(
      annotation = pipeline@annotation,
      genome_fa = pipeline@genome_fa,
      genome_bam = pipeline@genome_bam,
      outdir = pipeline@outdir,
      config = pipeline@config
    )
  }
  pipeline@novel_isoform_annotation <- novel_isoform_annotation
  annotation_to_fasta(
    isofornm_annotation = novel_isoform_annotation,
    genome_fa = pipeline@genome_fa,
    outfile = pipeline@transcriptome_assembly
  )
  return(pipeline)
})

setGeneric("read_realignment", function(pipeline) {
  standardGeneric("read_realignment")
})
setMethod("read_realignment", "FLAMES.Pipeline", function(pipeline) {
  if (pipeline@config$pipeline_parameters$oarfish_quantification) {
    minimap2_args <- c(
      "--eqx", "-N", "100", "-ax", "map-ont", "-y",
      "-t", pipeline@config$pipeline_parameters$threads,
      "--seed", pipeline@config$pipeline_parameters$seed
    )
  } else {
    minimap2_args <- c(
      "-ax", "map-ont", "-y", "-p", "0.9", "--end-bonus", "10", "-N", "3",
      "-t", pipeline@config$pipeline_parameters$threads,
      "--seed", pipeline@config$pipeline_parameters$seed 
    )
  }
  res <- lapply(
    seq_along(pipeline@fastq),
    function(i) {
      if (!is.null(names(pipeline@fastq))) {
        sample <- names(pipeline@fastq)[i]
      } else {
        sample <- pipeline@fastq[i]
      }
      message(sprintf("Realigning sample %s -> %s", sample, pipeline@transcriptome_bam[i]))
      minimap2_align(
        fq_in = pipeline@fastq[i],
        fa_file = pipeline@transcriptome_assembly,
        config = pipeline@config,
        outfile = pipeline@transcriptome_bam[i],
        minimap2_args = minimap2_args,
        minimap2 = pipeline@minimap2,
        samtools = pipeline@samtools,
        threads = pipeline@config$pipeline_parameters$threads
      )
    }
  )
  if (!is.null(names(pipeline@fastq))) {
    names(res) <- names(pipeline@fastq)
  }
  pipeline@metadata$read_realignment <- res
  return(pipeline)
})

setGeneric("transcript_quantification", function(pipeline) {
  standardGeneric("transcript_quantification")
})
setMethod("transcript_quantification", "FLAMES.Pipeline", function(pipeline, reference_only) {
  if ((!missing(reference_only) && reference_only) || is.na(pipeline@novel_isoform_annotation)) {
    annotation <- pipeline@annotation
  } else {
    annotation <- pipeline@novel_isoform_annotation
  }
  # TODO: refactor quantify_transcript to take file paths from pipeline
  pipeline@experiment <- quantify_transcript(
    annotation = annotation,
    outdir = pipeline@outdir,
    config = pipeline@config,
    pipeline = "bulk",
    samples = names(pipeline@fastq)
  )
  return(pipeline)
})

setGeneric("run_FLAMES", function(pipeline) {
  standardGeneric("run_FLAMES")
})

setMethod("run_FLAMES", "FLAMES.Pipeline", function(pipeline) {
  if (!prerun_check(pipeline, overwrite = FALSE)) {
    return(pipeline)
  }

  for (step in pipeline@steps) {
    # S4 objects are immutable
    # Need R6 for passing by reference
    pipeline <- run_step(pipeline, step)
    pipeline@completed_steps[step] <- TRUE
  }
  return(pipeline)
})

