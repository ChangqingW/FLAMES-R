setClass(
  "FLAMES.Pipeline",
  slots = list(
    # inputs
    config = "list",            # Configuration parameters
    outdir = "character",       # Output directory
    fastq = "character",        # Path to the FASTQ files
    annotation = "character",   # Path to the annotation file
    genome_fa = "character",    # Path to the genome FASTA file
    genome_mmi = "character",   # Path to the indexed genome file (minimap2)

    # outputs
    metadata = "list",          # Metadata for the pipeline
    bed = "character",         # Path to the BED file
    genome_bam = "character",     # Path to the genome BAM file
    transcriptome_bam = "character", # Path to the transcript BAM file
    novel_isoform_annotation = "character", # Path to the novel isoform GFF / GTF file
    transcriptome_assembly = "character", # Path to the transcriptome assembly file
    experiment = "SummarizedExperiment", # SummarizedExperiment object for quantification results

    # binaries
    minimap2 = "character",     # Path to the minimap2 binary
    samtools = "character",     # Path to the samtools binary

    # pipeline state
    steps = "logical",           # Steps to perform
    completed_steps = "logical", # Completed steps
    durations = "difftime",      # Durations of each step
    last_error = "list"
  ),
  prototype = list(
    config = list(),
    outdir = NA_character_,
    fastq = NA_character_,
    annotation = NA_character_,
    genome_fa = NA_character_,
    genome_mmi = NA_character_,
    metadata = list(),
    bed = NA_character_,
    genome_bam = NA_character_,
    transcriptome_bam = NA_character_,
    novel_isoform_annotation = NA_character_,
    transcriptome_assembly = NA_character_,
    # experiment = NA,

    minimap2 = NA_character_,
    samtools = NA_character_,
    steps = logical(),
    completed_steps = logical(),
    durations = structure(numeric(0), class = "difftime", units = "secs"),
    last_error = list()
  )
)

setClass(
  "FLAMES.SingleCellPipeline",
  contains = "FLAMES.Pipeline",
  slots = list(
    # inputs
    barcodes_file = "character",        # Path to the barcodes file
    expect_cell_number = "numeric",     # Expected number of cells
    # outputs
    demultiplexed_fastq = "character",  # path to demultiplexed FASTQ files
    deduped_fastq = "character"       # path to deduplicated FASTQ files
  ),
  prototype = list(
    barcodes_file = NA_character_,
    expect_cell_number = NA_real_,
    demultiplexed_fastq = NA_character_,
    deduped_fastq = NA_character_
  )
)

# Multi-sample pipeline is the same as SingleCellPipeline
# but input slots will be vectors of the same length
setClass(
  "FLAMES.MultiSampleSCPipeline",
  contains = "FLAMES.SingleCellPipeline",
  slots = list(
    experiments = "list" # List of SummarizedExperiment objects for each sample
  ),
  prototype = list(
    experiments = list()
  )
)

## Show methods
# helper method to display the pipeline class
setGeneric("display_pipeline_class", function(pipeline) {
  standardGeneric("display_pipeline_class")
})
#' @importFrom cli cli_alert cli_alert_warning cli_alert_success
setMethod("display_pipeline_class", "FLAMES.Pipeline", function(pipeline) {
  configured_steps <- pipeline@completed_steps[pipeline@steps]
  unfinished_steps <- names(which(!configured_steps))
  if (length(pipeline@last_error) > 0) {
    func <- cli::cli_alert_warning
  } else if (length(unfinished_steps) > 0) {
    func <- cli::cli_alert
  } else {
    func <- cli::cli_alert_success
  }
  func("A {.strong {class(pipeline)}} outputting to {.path {pipeline@outdir}}")
})

# Helper function to truncate paths
truncate_path <- function(path, max_chars = 50) {
  if (nchar(path) <= max_chars) {
    return(path)
  }
  paste0("...", substr(path, nchar(path) - max_chars + 4, nchar(path)))
}

# Helper function to format file sizes
format_file_size <- function(bytes) {
  units <- c("B", "KB", "MB", "GB", "TB")
  if (is.na(bytes) || bytes == 0) {
    return("0 B")
  }
  power <- min(floor(log(bytes, 1024)), length(units) - 1)
  sprintf("%.1f %s", bytes / (1024^power), units[power + 1])
}

#' @importFrom cli cli_h3 cli_alert_success cli_alert_warning cli_alert_info
display_inputs <- function(object, input_slots) {
  cli::cli_h3("Inputs")
  for (slot in input_slots) {
    paths <- slot(object, slot)
    if (length(paths) == 0 || (length(paths) == 1 && is.na(paths))) {
      if (slot == "barcodes_file" && !is.na(object@expect_cell_number)) {
        cli::cli_alert_info("{.field {slot}}: [not set] (set to expect {object@expect_cell_number} cells)")
      } else {
        cli::cli_alert_warning("{.field {slot}}: [not set]")
      }
    } else {
      if (length(paths) == 1) {
        if (file.exists(paths)) {
          cli::cli_alert_success("{.field {slot}}: {.path {truncate_path(paths)}}")
        } else {
          cli::cli_alert_warning("{.field {slot}}: {.path {truncate_path(paths)}} [missing]")
        }
      } else {
        if (all(!is.na(paths)) && all(file.exists(paths))) {
          cli::cli_alert_success("{.field {slot}}:")
        } else {
          cli::cli_alert_warning("{.field {slot}}:")
        }
        for (i in seq_along(paths)) {
          if (is.na(paths[i])) {
            if (slot == "genome_mmi") next
            if (slot == "barcodes_file" && !is.na(object@expect_cell_number[i])) {
              cli::cli_alert_info("  [not set] (set to expect {object@expect_cell_number[i]} cells)")
            } else {
              cli::cli_alert_warning("  [not set]")
            }
          } else if (file.exists(paths[i])) {
            cli::cli_alert_success("  {.path {truncate_path(paths[i])}}")
          } else {
            cli::cli_alert_warning("  {.path {truncate_path(paths[i])}} [missing]")
          }
        }
      }
    }
  }
}

#' @importFrom cli cli_h3 cli_alert_success cli_alert_info cli_text
display_outputs <- function(object, output_slots) {
  cli::cli_h3("Outputs")
  outdir <- object@outdir
  for (slot in output_slots) {
    paths <- slot(object, slot)
    if (length(paths) == 0 || (length(paths) == 1 && is.na(paths))) next

    display_path <- sapply(paths, function(path) {
      if (is.na(path)) {
        return(NA_character_)
      }
      tryCatch(
        {
          if (startsWith(
            normalizePath(path, mustWork = FALSE),
            normalizePath(outdir, mustWork = FALSE)
          )) {
            basename(path)
          } else {
            truncate_path(path)
          }
        },
        error = function(e) truncate_path(path)
      )
    })
    if (length(paths) == 1) {
      if (file.exists(paths[1])) {
        size <- format_file_size(file.info(paths[1])$size)
        cli::cli_alert_success("{.field {slot}}: {.path {display_path}} [{size}]")
      } else {
        cli::cli_alert_info("{.field {slot}}: {.path {display_path}}")
      }
    } else {
      cli::cli_text("{.field {slot}}:")
      for (i in seq_along(paths)) {
        if (is.na(paths[i])) {
          cli::cli_text("\t[not set]")
        } else if (file.exists(paths[i])) {
          size <- format_file_size(file.info(paths[i])$size)
          cli::cli_alert_success("  {.path {display_path[i]}} [{size}]")
        } else {
          cli::cli_alert_info("  {.path {display_path[i]}}")
        }
      }
    }
  }
}

# format difftime to human-friendly units
format_durations <- function(durations) {
  sapply(durations, function(d) {
    secs <- as.numeric(d, units = "secs")
    if (secs < 60) {
      sprintf("%.2f sec", secs)
    } else if (secs < 3600) {
      sprintf("%.2f min", secs / 60)
    } else if (secs < 86400) {
      sprintf("%.2f hr", secs / 3600)
    } else {
      sprintf("%.2f day", secs / 86400)
    }
  })
}

#' @importFrom cli cli_h3 cli_alert_success cli_alert_info cli_alert_danger cli_bullets
display_pipeline_steps <- function(object) {
  cli::cli_h3("Pipeline Steps")
  steps <- object@steps
  completed <- object@completed_steps
  durations <- object@durations

  for (step in names(steps)[steps]) { # configured steps only
    if (!is.null(completed[[step]]) && completed[[step]]) {
      duration <- durations[[step]]
      if (!is.null(duration)) {
        human_readable_duration <- format_durations(duration)
        cli::cli_alert_success("{.field {step}} (completed in {human_readable_duration})")
      } else {
        cli::cli_alert_success("{.field {step}} (completed)")
      }
    } else if (length(object@last_error) > 0 && step == object@last_error$step) {
      cli::cli_alert_danger("{.field {step}} (failed: {object@last_error$error})")
    } else {
      cli::cli_alert_info("{.field {step}} (pending)")
    }
  }
}

#' Show method for FLAMES.Pipeline
#'
#' Displays the pipeline in a pretty format
#' @param object An object of class `FLAMES.Pipeline`
#' @return None. Displays output to the console.
#' @examples
#' ppl <- example_pipeline()
#' show(ppl)
#' @importFrom methods show
#' @rdname show-FLAMESPipeline
#' @keywords internal
#' @export
setMethod("show", "FLAMES.Pipeline", function(object) {
  display_pipeline_class(object)
  display_inputs(object, c("fastq", "annotation", "genome_fa", "genome_mmi"))
  display_outputs(object, c(
    "genome_bam", "novel_isoform_annotation", "transcriptome_assembly", "transcriptome_bam"
  ))
  display_pipeline_steps(object)
})

#' @rdname show-FLAMESPipeline
#' @export
setMethod("show", "FLAMES.SingleCellPipeline", function(object) {
  display_pipeline_class(object)
  display_inputs(object, c("fastq", "annotation", "genome_fa", "genome_mmi", "barcodes_file"))
  display_outputs(object, c(
    "demultiplexed_fastq", "deduped_fastq",
    "genome_bam", "novel_isoform_annotation", "transcriptome_assembly", "transcriptome_bam"
  ))
  display_pipeline_steps(object)
})

#' @rdname show-FLAMESPipeline
#' @export
setMethod("show", "FLAMES.MultiSampleSCPipeline", function(object) {
  display_pipeline_class(object)
  display_inputs(object, c("fastq", "annotation", "genome_fa", "genome_mmi", "barcodes_file"))
  display_outputs(object, c(
    "novel_isoform_annotation", "transcriptome_assembly",
    "demultiplexed_fastq", "deduped_fastq",
    "genome_bam", "transcriptome_bam"
  ))
  display_pipeline_steps(object)
})

### getter and setters ###

#' Steps to perform in the pipeline
#'
#' @param object An object of class `FLAMES.Pipeline`
#' @return A named logical vector containing all possible steps
#' for the pipeline. The names of the vector are the step names,
#' and the values are logical indicating whether the step is
#' configured to be performed.
#' @examples
#' ppl <- example_pipeline()
#' steps(ppl)
#' @export
setGeneric("steps", function(object) standardGeneric("steps"))
#' @rdname steps
#' @export
setMethod("steps", "FLAMES.Pipeline", function(object) {
  object@steps
})

#' Set steps to perform in the pipeline
#'
#' @param object An object of class `FLAMES.Pipeline`
#' @param value A named logical vector containing all possible steps
#' for the pipeline. The names of the vector are the step names,
#' and the values are logical indicating whether the step is
#' configured to be performed.
#' @return An object of class `FLAMES.Pipeline` with the updated steps.
#' @examples
#' ppl <- example_pipeline()
#' steps(ppl) <- c(
#'   barcode_demultiplex = TRUE,
#'   genome_alignment = TRUE,
#'   gene_quantification = TRUE,
#'   isoform_identification = FALSE,
#'   read_realignment = FALSE,
#'   transcript_quantification = TRUE
#' )
#' ppl
#' # or partially change a step:
#' steps(ppl)["read_realignment"] <- TRUE
#' ppl
#' @export
setGeneric("steps<-", function(object, value) standardGeneric("steps<-"))
#' @rdname steps-set
#' @export
setMethod("steps<-", "FLAMES.Pipeline", function(object, value) {
  # validate the names
  if (any(is.null(names(value))) || any(is.na(names(value)))) {
    stop("Steps must be a named logical vector.")
  }
  if (!all(names(value) %in% names(object@steps))) {
    stop(sprintf(
      "Invalid step names. Expected: %s, but got: %s",
      paste(names(object@steps), collapse = ", "),
      paste(names(value)[!names(value) %in% names(object@steps)], collapse = ", ")
    ))
  }
  object@steps[names(value)] <- value
  object
})
