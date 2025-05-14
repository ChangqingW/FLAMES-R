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
  } else if (!all(utils::file_test("-f", fastq))) {
    stop("fastq must be a valid path to a folder or a FASTQ file")
  }
  if (is.null(names(fastq))) {
    names(fastq) <- gsub("\\.(fastq|fq)(\\.gz)?$", "", basename(fastq))
  }
  names(fastq) <- make.unique(names(fastq), sep = "_")

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
  pipeline@transcriptome_bam <- file.path(outdir, paste0(names(fastq), "_", "realign2transcript.bam"))
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
  ##

  ## pipeline state
  pipeline@steps <- steps
  pipeline@completed_steps <- setNames(
    rep(FALSE, length(steps)), names(steps)
  )
  # TODO: add resume option
  # validate if e.g. genome_bam exists, skip genome alignment step

  return(pipeline)
}

setGeneric("prerun_check", function(pipeline, overwrite = FALSE) {
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
  start_time <- Sys.time()
  message(sprintf("Running step: %s", step))
  pipeline <- switch(step,
    barcode_demultiplex = barcode_demultiplex(pipeline),
    genome_alignment = genome_alignment(pipeline),
    gene_quantification = gene_quantification(pipeline),
    isoform_identification = isoform_identification(pipeline),
    read_realignment = read_realignment(pipeline),
    transcript_quantification = transcript_quantification(pipeline),
    stop(sprintf("Unknown step: %s", step))
  )
  end_time <- Sys.time()
  pipeline@completed_steps[step] <- TRUE
  pipeline@durations[step] <- difftime(end_time, start_time, units = "secs")
  return(pipeline)
})

# individual steps as methods
setGeneric("barcode_demultiplex", function(pipeline) {
  standardGeneric("barcode_demultiplex")
})
setMethod("barcode_demultiplex", "FLAMES.Pipeline", function(pipeline) {
  stop("Barcode demultiplexing is only implemented for single cell pipelines (SingleCellPipeline() and MultiSampleSCPipeline())")
})
setGeneric("gene_quantification", function(pipeline) {
  standardGeneric("gene_quantification")
})
setMethod("gene_quantification", "FLAMES.Pipeline", function(pipeline) {
  # todo: implement gene quantification for bulk pipelines
  stop("Gene quantification is not implemented for bulk pipelines yet")
})
setGeneric("genome_alignment_raw", function(pipeline, fastqs) {
  standardGeneric("genome_alignment_raw")
})
setMethod("genome_alignment_raw", "FLAMES.Pipeline", function(pipeline, fastqs) {
  minimap2_args <- c(
    "-ax", "splice", "-k14", "--secondary=no", # "-y",
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
    seq_along(fastqs),
    function(i) {
      if (!is.null(names(fastqs))) {
        sample <- names(fastqs)[i]
      } else {
        sample <- fastqs[i]
      }
      message(sprintf("Aligning sample %s -> %s", sample, pipeline@genome_bam[i]))
      minimap2_align(
        fq_in = fastqs[i],
        fa_file = pipeline@genome_fa,
        config = pipeline@config,
        outfile = pipeline@genome_bam[i],
        minimap2_args = minimap2_args,
        sort_by = "coordinates",
        minimap2 = pipeline@minimap2,
        samtools = pipeline@samtools,
        threads = pipeline@config$pipeline_parameters$threads
      )
    }
  )
  if (!is.null(names(fastqs))) {
    names(res) <- names(fastqs)
  }
  unlink(bed_file)
  pipeline@metadata$genome_alignment <- res
  return(pipeline)
})
setGeneric("genome_alignment", function(pipeline) {
  standardGeneric("genome_alignment")
})
setMethod("genome_alignment", "FLAMES.Pipeline", function(pipeline) {
  genome_alignment_raw(
    pipeline = pipeline,
    fastqs = pipeline@fastq
  )
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
    isoform_annotation = novel_isoform_annotation,
    genome_fa = pipeline@genome_fa,
    outfile = pipeline@transcriptome_assembly
  )
  return(pipeline)
})

setGeneric("read_realignment_raw", function(pipeline, include_tags, sort_by, fastqs) {
  standardGeneric("read_realignment_raw")
})
setMethod(
  "read_realignment_raw", "FLAMES.Pipeline",
  function(pipeline, include_tags = FALSE, sort_by, fastqs) {
    if (pipeline@config$pipeline_parameters$oarfish_quantification) {
      minimap2_args <- c(
        "--eqx", "-N", "100", "-ax", "map-ont",
        "-t", pipeline@config$pipeline_parameters$threads,
        "--seed", pipeline@config$pipeline_parameters$seed
      )
    } else {
      minimap2_args <- c(
        "-ax", "map-ont", "-p", "0.9", "--end-bonus", "10", "-N", "3",
        "-t", pipeline@config$pipeline_parameters$threads,
        "--seed", pipeline@config$pipeline_parameters$seed
      )
    }
    if (include_tags) {
      minimap2_args <- base::append(minimap2_args, "-y")
    }

    res <- lapply(
      seq_along(fastqs),
      function(i) {
        if (!is.null(names(fastqs))) {
          sample <- names(fastqs)[i]
        } else {
          sample <- fastqs[i]
        }
        message(sprintf("Realigning sample %s -> %s", sample, pipeline@transcriptome_bam[i]))
        minimap2_align(
          fq_in = fastqs[i],
          fa_file = pipeline@transcriptome_assembly,
          config = pipeline@config,
          outfile = pipeline@transcriptome_bam[i],
          minimap2_args = minimap2_args,
          sort_by = sort_by,
          minimap2 = pipeline@minimap2,
          samtools = pipeline@samtools,
          threads = pipeline@config$pipeline_parameters$threads
        )
      }
    )
    if (!is.null(names(fastqs))) {
      names(res) <- names(fastqs)
    }
    pipeline@metadata$read_realignment <- res
    return(pipeline)
  }
)

setGeneric("read_realignment", function(pipeline, include_tags) {
  standardGeneric("read_realignment")
})
setMethod("read_realignment", "FLAMES.Pipeline", function(pipeline, include_tags = FALSE) {
  sort_by <- ifelse(
    pipeline@config$pipeline_parameters$oarfish_quantification,
    "none",
    "coordinates"
  )
  read_realignment_raw(
    pipeline = pipeline,
    include_tags = include_tags,
    sort_by = sort_by,
    fastqs = pipeline@fastq
  )
})
setGeneric("transcript_quantification", function(pipeline, reference_only) {
  standardGeneric("transcript_quantification")
})
setMethod("transcript_quantification", "FLAMES.Pipeline", function(pipeline, reference_only) {
  if ((!missing(reference_only) && reference_only) || is.na(pipeline@novel_isoform_annotation)) {
    annotation <- pipeline@annotation
  } else {
    annotation <- pipeline@novel_isoform_annotation
  }
  pipeline_class <- switch(
    class(pipeline),
    "FLAMES.Pipeline" = "bulk",
    "FLAMES.SingleCellPipeline" = "sc_single_sample",
    "FLAMES.MultiSampleSCPipeline" = "sc_multi_sample"
  )
  # TODO: refactor quantify_transcript to take file paths from pipeline
  x <- quantify_transcript(
    annotation = annotation,
    outdir = pipeline@outdir,
    config = pipeline@config,
    pipeline = pipeline_class,
    samples = names(pipeline@fastq)
  )
  if (is.list(x) & pipeline_class == "sc_multi_sample") {
    pipeline@experiments <- x
  } else {
    pipeline@experiment <- x
  }
  return(pipeline)
})

setGeneric("run_FLAMES", function(pipeline) {
  standardGeneric("run_FLAMES")
})
setMethod("run_FLAMES", "FLAMES.Pipeline", function(pipeline) {
  if (!prerun_check(pipeline, overwrite = FALSE)) {
    return(pipeline)
  }

  for (step in names(which(pipeline@steps))) {
    # S4 objects are immutable
    # Need R6 for passing by reference
    pipeline <- tryCatch(
      run_step(pipeline, step),
      error = function(e) {
        warning(sprintf("Error in step %s: %s, pipeline stopped.", step, e$message))
        pipeline@last_error <- list(
          step = step,
          error = e,
          traceback = capture.output(traceback())
        )
        return(pipeline)
      }
    )
    pipeline@completed_steps[step] <- TRUE
  }
  return(pipeline)
})

setGeneric("resume_FLAMES", function(pipeline) {
  standardGeneric("resume_FLAMES")
})
setMethod("resume_FLAMES", "FLAMES.Pipeline", function(pipeline) {
  configured_steps <- pipeline@completed_steps[pipeline@steps]
  unfinished_steps <- names(which(!configured_steps))
  if (length(unfinished_steps) == 0) {
    message("All steps have already been completed.")
    return(pipeline)
  } else {
    message("Resuming pipeline from step: ", unfinished_steps[1])
    for (step in unfinished_steps) {
      pipeline <- tryCatch(
        run_step(pipeline, step),
        error = function(e) {
          warning(sprintf("Error in step %s: %s, pipeline stopped.", step, e$message))
          pipeline@last_error <- list(
            step = step,
            error = e,
            traceback = capture.output(traceback())
          )
          return(pipeline)
        }
      )
      pipeline@completed_steps[step] <- TRUE
    }
  }
})
