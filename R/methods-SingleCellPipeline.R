#' Pipeline for Single Cell Data
#'
#' @md
#'
#' @description
#' Semi-supervised isofrom detection and annotation for long read data. This variant is
#' meant for single sample scRNA-seq data. Specific parameters can be configured in
#' the config file (see \code{\link{create_config}}), input files are specified via
#' arguments.
#'
#' @details
#' By default the pipeline starts with demultiplexing the input fastq data. If the
#' cell barcodes are known apriori (e.g. via coupled short-read sequencing), the
#' \code{barcodes_file} argument can be used to specify a file containing the cell
#' barcodes, and a modified Rcpp version of \code{flexiplex} will be used; otherwise,
#' \code{expect_cell_number} need to be provided, and \code{BLAZE} will be used to
#' generate the cell barcodes. The pipeline then aligns the reads to the genome using
#' \code{minimap2}. The alignment is then used for isoform detection (either using
#' \code{FLAMES} or \code{bambu}, can be configured). The reads are then realigned
#' to the detected isoforms. Finally, a transcript count matrix is generated (either
#' using \code{FLAMES}'s simplistic counting or \code{oarfish}'s Expectation
#' Maximization algorithm, can be configured). The results can be accssed with
#' \code{experiment(pipeline)}. If the pipeline errored out / new steps were configured,
#' it can be resumed by calling \code{resume_FLAMES(pipeline)}
#'
#' @inheritParams BulkPipeline
#' @param barcodes_file The file with expected cell barcodes, with each barcode on a new line.
#' @param expect_cell_number The expected number of cells in the sample. This is used if
#'   \code{barcodes_file} is not provided. See \code{BLAZE} for more details.
#'
#' @return A \code{FLAMES.SingleCellPipeline} object. The pipeline can be run using
#' \code{run_FLAMES(pipeline)}. The results can be accessed with \code{experiment(pipeline)}.
#' The pipeline also outputs a number of output files into the given \code{outdir} directory.
#' Some of these output files include:
#' \describe{
#'  \item{matched_reads.fastq}{ - fastq file with reads demultiplexed}
#'  \item{align2genome.bam}{ - sorted BAM file with reads aligned to genome}
#'  \item{matched_reads_dedup.fastq}{ - demultiplexed and UMI-deduplicated fastq file}
#'  \item{transcript_assembly.fa}{ - transcript sequence from the isoforms}
#'  \item{isoform_annotated.filtered.gff3}{ - isoforms in gff3 format (also contained in the SingleCellExperiment)}
#'  \item{realign2transcript.bam}{ - sorted realigned BAM file using the transcript_assembly.fa as reference}
#' }
#'
#' @seealso
#' \code{\link{create_config}} for creating a configuration file,
#' \code{\link{BulkPipeline}} for bulk long data,
#' \code{\link{MultiSampleSCPipeline}} for multi sample single cell pipelines.
#'
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' bc_allow <- file.path(outdir, "bc_allow.tsv")
#' genome_fa <- file.path(outdir, "rps24.fa")
#' R.utils::gunzip(
#'   filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
#'   destname = bc_allow, remove = FALSE
#' )
#' R.utils::gunzip(
#'   filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
#'   destname = genome_fa, remove = FALSE
#' )
#' ppl <- SingleCellPipeline(
#'   config_file = create_config(outdir, gene_quantification = FALSE),
#'   outdir = outdir,
#'   fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
#'   annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
#'   genome_fa = genome_fa,
#'   barcodes_file = bc_allow
#' )
#' ppl <- run_FLAMES(ppl)
#' experiment(ppl)
#'
#' @export
SingleCellPipeline <- function(
  config_file, outdir, fastq, annotation, genome_fa, genome_mmi,
  minimap2, samtools, barcodes_file, expect_cell_number, controllers
) {
  pipeline <- new("FLAMES.SingleCellPipeline")
  config <- check_arguments(annotation, fastq, genome_bam = NULL, outdir, genome_fa, config_file)$config

  if (!dir.exists(outdir)) {
    dir.create(outdir)
    message(sprintf("Output directory (%s) did not exist, created.", outdir))
  }

  steps <- c(
    "barcode_demultiplex", "genome_alignment", "gene_quantification",
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
  if (!missing(genome_mmi) && is.character(genome_mmi)) {
    pipeline@genome_mmi <- genome_mmi
  }
  if (!missing(barcodes_file) && is.character(barcodes_file)) {
    pipeline@barcodes_file <- barcodes_file
  } else if (!missing(expect_cell_number) && is.numeric(expect_cell_number)) {
    pipeline@expect_cell_number <- expect_cell_number
  } else if (steps["barcode_demultiplex"]) {
    stop("Either barcodes_file or expect_cell_number must be provided.")
  }

  ## outputs
  # metadata
  pipeline@bed <- file.path(outdir, "reference.bed")
  pipeline@genome_bam <- file.path(outdir, "align2genome.bam")
  pipeline@transcriptome_bam <- file.path(outdir, "realign2transcript.bam")
  pipeline@transcriptome_assembly <- file.path(outdir, "transcript_assembly.fa")
  pipeline@demultiplexed_fastq <- file.path(outdir, "matched_reads.fastq.gz")
  pipeline@deduped_fastq <- file.path(outdir, "matched_reads_dedup.fastq.gz")

  ## binaries
  if (missing(minimap2) || !is.character(minimap2)) {
    minimap2 <- find_bin("minimap2")
    if (is.na(minimap2)) {
      stop("minimap2 not found, please make sure it is installed and provide its path as the minimap2 argument")
    }
  }
  pipeline@minimap2 <- minimap2
  if (missing(samtools) || !is.character(samtools)) {
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

  if (!missing(controllers)) {
    pipeline@controllers <- normalize_controllers(controllers, names(steps))
  }

  if (pipeline@config$pipeline_parameters$multithread_isoform_identification &
    !pipeline@config$pipeline_parameters$bambu_isoform_identification) {
    warning("The multithreaded isoform identification implmentation is currently unstable and may throw errors. Set `multithread_isoform_identification` to `FALSE` in the config file or with `pipeline@config$pipeline_parameters$multithread_isoform_identification <- FALSE` to fall back to the single-threaded implementation. Report to https://github.com/mritchielab/FLAMES/issues if you encounter any problems.")
  }
  return(pipeline)
}

#' Example pipelins
#'
#' @description
#' Provides example pipelines for bulk, single cell and multi-sample single cell.
#'
#' @param type The type of pipeline to create. Options are "SingleCellPipeline",
#'   "BulkPipeline", and "MultiSampleSCPipeline".
#' @param outdir (Optional) The output directory where the example pipeline will
#'   be created. If not provided, a temporary directory will be created.
#'
#' @return A pipeline object of the specified type.
#'
#' @seealso
#' \code{\link{SingleCellPipeline}} for creating the single cell pipeline,
#' \code{\link{BulkPipeline}} for bulk long data,
#' \code{\link{MultiSampleSCPipeline}} for multi sample single cell pipelines.
#'
#' @examples
#' example_pipeline("SingleCellPipeline")
#'
#' @importFrom R.utils gunzip
#' @importFrom ShortRead readFastq writeFastq
#'
#' @export
example_pipeline <- function(type = "SingleCellPipeline", outdir) {
  if (missing(outdir)) {
    outdir <- tempfile()
    dir.create(outdir)
  }
  switch(type,
    "SingleCellPipeline" = {
      bc_allow <- file.path(outdir, "bc_allow.tsv")
      genome_fa <- file.path(outdir, "rps24.fa")
      R.utils::gunzip(
        filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
        destname = bc_allow, remove = FALSE
      )
      R.utils::gunzip(
        filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
        destname = genome_fa, remove = FALSE
      )
      SingleCellPipeline(
        config_file = create_config(outdir, gene_quantification = FALSE),
        outdir = outdir,
        fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
        annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
        genome_fa = genome_fa,
        barcodes_file = bc_allow
      )
    },
    "BulkPipeline" = {
      reads <- ShortRead::readFastq(
        system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES")
      )
      dir.create(file.path(outdir, "fastq"))
      genome_fa <- file.path(outdir, "rps24.fa")
      R.utils::gunzip(
        filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
        destname = genome_fa, remove = FALSE
      )
      ShortRead::writeFastq(reads[1:100],
        file.path(outdir, "fastq/sample1.fq.gz"), mode = "w", full = FALSE)
      reads <- reads[-(1:100)]
      ShortRead::writeFastq(reads[1:100],
        file.path(outdir, "fastq/sample2.fq.gz"), mode = "w", full = FALSE)
      reads <- reads[-(1:100)]
      ShortRead::writeFastq(reads,
        file.path(outdir, "fastq/sample3.fq.gz"), mode = "w", full = FALSE)
      BulkPipeline(
        fastq = c(
          "sample1" = file.path(outdir, "fastq", "sample1.fq.gz"),
          "sample2" = file.path(outdir, "fastq", "sample2.fq.gz"),
          "sample3" = file.path(outdir, "fastq", "sample3.fq.gz")
        ),
        annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
        genome_fa = genome_fa,
        config_file = create_config(outdir, type = "sc_3end", threads = 1, no_flank = TRUE),
        outdir = outdir
      )
    },
    "MultiSampleSCPipeline" = {
      reads <- ShortRead::readFastq(
        system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES")
      )
      dir.create(file.path(outdir, "fastq"))
      bc_allow <- file.path(outdir, "bc_allow.tsv")
      genome_fa <- file.path(outdir, "rps24.fa")
      R.utils::gunzip(
        filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
        destname = bc_allow, remove = FALSE
      )
      R.utils::gunzip(
        filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
        destname = genome_fa, remove = FALSE
      )
      ShortRead::writeFastq(reads[1:100],
        file.path(outdir, "fastq/sample1.fq.gz"), mode = "w", full = FALSE)
      reads <- reads[-(1:100)]
      ShortRead::writeFastq(reads[1:100],
        file.path(outdir, "fastq/sample2.fq.gz"), mode = "w", full = FALSE)
      reads <- reads[-(1:100)]
      ShortRead::writeFastq(reads,
        file.path(outdir, "fastq/sample3.fq.gz"), mode = "w", full = FALSE)
      MultiSampleSCPipeline(
        config_file = create_config(outdir, type = "sc_3end", threads = 1, no_flank = TRUE),
        outdir = outdir,
        fastq = c("sampleA" = file.path(outdir, "fastq"),
          "sample1" = file.path(outdir, "fastq", "sample1.fq.gz"),
          "sample2" = file.path(outdir, "fastq", "sample2.fq.gz"),
          "sample3" = file.path(outdir, "fastq", "sample3.fq.gz")),
        annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
        genome_fa = genome_fa,
        barcodes_file = rep(bc_allow, 4),
        controllers = crew::crew_controller_local()
      )
    }
  )
}

setMethod("barcode_demultiplex", "FLAMES.SingleCellPipeline", function(pipeline) {
  if (any(is.na(pipeline@barcodes_file)) && any(is.na(pipeline@expect_cell_number))) {
    stop("Either barcodes_file or expect_cell_number must be provided.")
  }

  using_controllers <- FALSE
  if ("barcode_demultiplex" %in% names(pipeline@controllers)) {
    using_controllers <- TRUE
    controller <- pipeline@controllers[["barcode_demultiplex"]]
  } else if ("default" %in% names(pipeline@controllers)) {
    using_controllers <- TRUE
    controller <- pipeline@controllers[["default"]]
  }

  if (any(is.na(pipeline@barcodes_file))) {
    message("No barcodes file provided, running BLAZE to generate barcode list from long reads...")
    if (using_controllers) {
      controller$start()
      controller$push(
        command = blaze(
          expect_cells = expect_cell_number,
          fq_in = fastq,
          "output-prefix" = file.path(outdir, ""),
          "output-fastq" = demultiplexed_fastq,
          "threads" = config$pipeline_parameters$threads,
          "max-edit-distance" = config$barcode_parameters$max_bc_editdistance,
          "overwrite" = TRUE
        ),
        data = list(
          expect_cell_number = pipeline@expect_cell_number,
          fastq = pipeline@fastq,
          outdir = pipeline@outdir,
          demultiplexed_fastq = pipeline@demultiplexed_fastq,
          config = pipeline@config,
          blaze = FLAMES:::blaze
        )
      )
      controller$wait(mode = "all")
      task <- controller$pop(error = "stop")
      # res <- task$result[[1]]
      controller$terminate()
    } else {
      blaze(
        expect_cells = pipeline@expect_cell_number,
        fq_in = pipeline@fastq,
        "output-prefix" = file.path(pipeline@outdir, ""),
        "output-fastq" = pipeline@demultiplexed_fastq,
        "threads" = pipeline@config$pipeline_parameters$threads,
        "max-edit-distance" = pipeline@config$barcode_parameters$max_bc_editdistance,
        "overwrite" = TRUE
      )
    }
  } else {
    if (using_controllers) {
      controller$start()
      controller$push(
        command = find_barcode(
          fastq = fastq,
          barcodes_file = barcodes_file,
          stats_out = file.path(outdir, "matched_barcode_stat"),
          reads_out = demultiplexed_fastq,
          pattern = setNames(
            as.character(config$barcode_parameters$pattern),
            names(config$barcode_parameters$pattern)
          ),
          TSO_seq = config$barcode_parameters$TSO_seq,
          TSO_prime = config$barcode_parameters$TSO_prime,
          cutadapt_minimum_length = config$barcode_parameters$cutadapt_minimum_length,
          full_length_only = config$barcode_parameters$full_length_only,
          max_bc_editdistance = config$barcode_parameters$max_bc_editdistance,
          max_flank_editdistance = config$barcode_parameters$max_flank_editdistance,
          strand = config$barcode_parameters$strand,
          threads = config$pipeline_parameters$threads
        ),
        data = list(
          fastq = pipeline@fastq,
          barcodes_file = pipeline@barcodes_file,
          outdir = pipeline@outdir,
          demultiplexed_fastq = pipeline@demultiplexed_fastq,
          config = pipeline@config,
          find_barcode = FLAMES:::find_barcode
        )
      )
      controller$wait(mode = "all")
      task <- controller$pop(error = "stop")
      res <- task$result[[1]]
      controller$terminate()
    } else {
      res <- find_barcode(
        fastq = pipeline@fastq,
        barcodes_file = pipeline@barcodes_file,
        stats_out = file.path(pipeline@outdir, "matched_barcode_stat"),
        reads_out = pipeline@demultiplexed_fastq,
        pattern = setNames(
          as.character(pipeline@config$barcode_parameters$pattern),
          names(pipeline@config$barcode_parameters$pattern)
        ),
        TSO_seq = pipeline@config$barcode_parameters$TSO_seq,
        TSO_prime = pipeline@config$barcode_parameters$TSO_prime,
        cutadapt_minimum_length = pipeline@config$barcode_parameters$cutadapt_minimum_length,
        full_length_only = pipeline@config$barcode_parameters$full_length_only,
        max_bc_editdistance = pipeline@config$barcode_parameters$max_bc_editdistance,
        max_flank_editdistance = pipeline@config$barcode_parameters$max_flank_editdistance,
        strand = pipeline@config$barcode_parameters$strand,
        threads = pipeline@config$pipeline_parameters$threads
      )
    }
    pipeline@metadata$barcode_demultiplex <- res
  }
  return(pipeline)
})

setMethod("genome_alignment", "FLAMES.SingleCellPipeline", function(pipeline) {
  infq <- pipeline@fastq
  if (all(file.exists(pipeline@demultiplexed_fastq))) {
    infq <- pipeline@demultiplexed_fastq
  }
  genome_alignment_raw(
    pipeline = pipeline,
    fastqs = infq
  )
})

setMethod("gene_quantification", "FLAMES.SingleCellPipeline", function(pipeline) {
  infq <- pipeline@fastq
  if (all(file.exists(pipeline@demultiplexed_fastq))) {
    infq <- pipeline@demultiplexed_fastq
  }
  pipeline_class <- switch(
    class(pipeline),
    "FLAMES.SingleCellPipeline" = "sc_single_sample",
    "FLAMES.MultiSampleSCPipeline" = "sc_multi_sample"
  )
  quantify_gene(
    annotation = pipeline@annotation,
    outdir = pipeline@outdir,
    pipeline = pipeline_class,
    infq = infq,
    in_bam = pipeline@genome_bam,
    out_fastq = pipeline@deduped_fastq,
    n_process = pipeline@config$pipeline_parameters$threads,
    saturation_curve = TRUE,
    sample_names = names(pipeline@fastq),
    random_seed = pipeline@config$pipeline_parameters$seed
  )
  return(pipeline)
})
setMethod("read_realignment", "FLAMES.SingleCellPipeline", function(pipeline, include_tags = FALSE) {
  # work out which fastq to use
  for (fastq in list(
    pipeline@fastq,
    pipeline@demultiplexed_fastq,
    pipeline@deduped_fastq
  )) {
    message(sprintf("Checking for fastq file(s) %s", paste(fastq, collapse = ", ")))
    if (all(file.exists(fastq))) {
      reads_to_realign <- fastq
      message("\tfiles found")
    } else {
      message("\tfiles not found")
    }
  }
  for_oarfish <- pipeline@config$pipeline_parameters$oarfish_quantification
  if (for_oarfish) {
    if (!all(file.exists(pipeline@deduped_fastq))) {
      warning("Oarfish does not support UMI deduplication, you should deduplicate reads before running Oarfish")
    }
    if (!missing(include_tags) && !include_tags) {
      warning("Oarfish need CB tags for single-cell quantification but include_tags is set to FALSE")
    } else {
      include_tags <- TRUE
    }
    sort_by <- "CB"
  } else {
    sort_by <- "coordinates"
  }
  read_realignment_raw(
    pipeline = pipeline,
    include_tags = include_tags,
    sort_by = sort_by,
    fastqs = reads_to_realign
  )
})

#' Plot Cell Barcode demultiplex statistics
#'
#' @description produce a barplot of cell barcode demultiplex statistics
#' @param pipeline A \code{FLAMES.SingleCellPipeline} object
#' @return a list of ggplot objects:
#' \itemize{
#' \item reads_count_plot: stacked barplot of: demultiplexed reads
#' \item knee_plot: knee plot of UMI counts before TSO trimming
#' \item flank_editdistance_plot: flanking sequence (adaptor) edit-distance plot
#' \item barcode_editdistance_plot: barcode edit-distance plot
#' \item cutadapt_plot: if TSO trimming is performed, number of reads kept by cutadapt
#' }
#' @examples
#' pipeline <- example_pipeline("MultiSampleSCPipeline") |>
#'   run_step("barcode_demultiplex")
#' plot_demultiplex(pipeline)
#' @rdname plot_demultiplex
#' @export
setGeneric("plot_demultiplex", function(pipeline) {
  standardGeneric("plot_demultiplex")
})
#' @rdname plot_demultiplex
#' @export
setMethod("plot_demultiplex", "FLAMES.SingleCellPipeline", function(pipeline) {
  if (is.null(pipeline@metadata$barcode_demultiplex)) {
    stop("No barcode demultiplexing results found, have you run the demultiplex step?")
  }
  plot_demultiplex_raw(pipeline@metadata$barcode_demultiplex)
})

#' Pipeline for Single Cell Data (deprecated)
#'
#' @description This function is deprecated. Please use [SingleCellPipeline()] instead.
#'
#' @param annotation The file path to the annotation file in GFF3 format
#' @param fastq The file path to input fastq file
#' @param outdir The path to directory to store all output files.
#' @param genome_fa The file path to genome fasta file.
#' @param minimap2 Path to minimap2, optional.
#' @param config_file File path to the JSON configuration file.
#' @param barcodes_file The file with expected cell barcodes, with each barcode on a new line.
#' @param expect_cell_number The expected number of cells in the sample. This is used if
#'   \code{barcodes_file} is not provided. See \code{BLAZE} for more details.
#'
#' @return A \code{SingleCellPipeline} object containing the transcript counts.
#'
#' @seealso
#' \code{\link{SingleCellPipeline}} for the new pipeline interface,
#' \code{\link{BulkPipeline}} for bulk long data,
#' \code{\link{MultiSampleSCPipeline}} for multi sample single cell pipelines.
#'
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' bc_allow <- file.path(outdir, "bc_allow.tsv")
#' genome_fa <- file.path(outdir, "rps24.fa")
#' R.utils::gunzip(
#'   filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
#'   destname = bc_allow, remove = FALSE
#' )
#' R.utils::gunzip(
#'   filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
#'   destname = genome_fa, remove = FALSE
#' )
#' sce <- FLAMES::sc_long_pipeline(
#'   genome_fa = genome_fa,
#'   fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
#'   annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
#'   outdir = outdir,
#'   barcodes_file = bc_allow
#' )
#' @export
sc_long_pipeline <- function(
    annotation, fastq, outdir, genome_fa, minimap2 = NULL,
    barcodes_file = NULL, expect_cell_number = NULL, config_file = NULL) {
  pipeline <- SingleCellPipeline(
    config_file = config_file,
    outdir = outdir,
    fastq = fastq,
    annotation = annotation,
    genome_fa = genome_fa,
    minimap2 = minimap2,
    barcodes_file = barcodes_file,
    expect_cell_number = expect_cell_number
  )
  pipeline <- run_FLAMES(pipeline)
  saveRDS(pipeline, file.path(outdir, "pipeline.rds"))
  message("Pipeline saved to ", file.path(outdir, "pipeline.rds"))
  if (length(pipeline@last_error) == 0) {
    return(experiment(pipeline))
  } else {
    warning("Returning pipeline object instead of experiment due to errors.")
    message("You can resume the pipeline after resolving the errors with resume_FLAMES(pipeline)")
    return(pipeline)
  }
}

#' @importFrom utils read.csv
#' @importFrom dplyr summarise_at group_by
#' @importFrom GenomicRanges GRangesList GRanges
generate_sc_sce <- function(out_files, create_function) {

  transcript_count <- read.csv(out_files$counts, stringsAsFactors = FALSE)
  if ("fsm_annotation" %in% names(out_files)) {
    isoform_FSM_annotation <- read.csv(out_files$fsm_annotation, stringsAsFactors = FALSE)
  } else {
    isoform_FSM_annotation <- read.csv(file.path(out_files$outdir, "isoform_FSM_annotation.csv"), stringsAsFactors = FALSE)
  }
  transcript_count <- transcript_count[transcript_count$transcript_id %in% isoform_FSM_annotation$transcript_id, ]

  isoform_FSM_annotation <- isoform_FSM_annotation[match(transcript_count$transcript_id, isoform_FSM_annotation$transcript_id), ]
  transcript_count$FSM_match <- isoform_FSM_annotation$FSM_match
  if (!all(transcript_count$transcript_id %in% isoform_FSM_annotation$transcript_id)) {
    message("Some transcript_ids are not recorded in isoform_FSM_annotation.csv")
    transcript_count$FSM_match[is.na(transcript_count$FSM_match)] <- transcript_count$transcript_id[is.na(transcript_count$FSM_match)]
  }

  cell_bcs <- colnames(transcript_count)[!(colnames(transcript_count) %in% c("transcript_id", "gene_id", "FSM_match"))]
  tr_anno <- transcript_count[, c("transcript_id", "gene_id", "FSM_match")]

  # sum transcript (FSM) counts
  mer_tmp <- transcript_count |>
    dplyr::group_by(FSM_match) |> # nolint: object_usage_linter.
    dplyr::summarise_at(cell_bcs, sum)

  # Create long read SCE
  tr_anno <- tr_anno[match(mer_tmp$FSM_match, tr_anno$FSM_match), ]
  tr_sce <- create_function(
    assays = list(counts = as.matrix(mer_tmp[, -1])),
    metadata = list("OutputFiles" = out_files),
  )

  rowData(tr_sce) <- DataFrame(tr_anno)
  rownames(tr_sce) <- tr_anno$FSM_match
  return(tr_sce)
}

generate_sc_singlecell <- function(out_files) {
  return(generate_sc_sce(out_files = out_files, create_function = SingleCellExperiment::SingleCellExperiment))
}

generate_bulk_summarized <- function(out_files) {
  return(generate_sc_sce(out_files = out_files, create_function = SummarizedExperiment::SummarizedExperiment))
}

#' Create \code{SingleCellExperiment} object from \code{FLAMES} output folder
#' @param outdir The folder containing \code{FLAMES} output files
#' @param annotation the annotation file that was used to produce the output files
#' @param quantification (Optional)  the quantification method used to generate the 
#' output files (either "FLAMES" or "Oarfish".). If not specified, the function will 
#' attempt to determine the quantification method.
#' @return a list of \code{SingleCellExperiment} objects if multiple transcript matrices were
#' found in the output folder, or a \code{SingleCellExperiment} object if only one were found
#' @export
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' bc_allow <- file.path(outdir, "bc_allow.tsv")
#' genome_fa <- file.path(outdir, "rps24.fa")
#' R.utils::gunzip(
#'   filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
#'   destname = bc_allow, remove = FALSE
#' )
#' R.utils::gunzip(
#'   filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
#'   destname = genome_fa, remove = FALSE
#' )
#' annotation <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
#'
#' sce <- sc_long_pipeline(
#'   genome_fa = genome_fa,
#'   fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
#'   annotation = annotation,
#'   outdir = outdir,
#'   barcodes_file = bc_allow,
#'   config_file = create_config(outdir, oarfish_quantification = FALSE)
#' )
#' sce_2 <- create_sce_from_dir(outdir, annotation)
create_sce_from_dir <- function(outdir, annotation, quantification = "FLAMES") {

  samples <- list.files(outdir, pattern = "_?transcript_count.csv.gz$", full.names = TRUE, all.files = TRUE)
  samples_oarfish <- list.files(outdir, pattern = "\\.count\\.mtx$", full.names = TRUE, all.files = TRUE) |>
    stringr::str_remove("\\.count\\.mtx$")

  if (length(samples) > 0 && quantification == "FLAMES") {
    sce_list <- lapply(samples, \(x) {
      out_files <- list(
        counts = x,
        outdir = outdir,
        transcript_assembly = file.path(outdir, "transcript_assembly.fa"),
        annotation = annotation
      )
      generate_sc_singlecell(out_files) |>
        addRowRanges(annotation, outdir)
    }) |>
      setNames(stringr::str_remove(samples, "_?transcript_count.csv.gz$"))

  } else if (length(samples_oarfish) > 0 && (missing(quantification) || quantification == "Oarfish")) {
    sce_list <- lapply(samples_oarfish, \(x) {
      parse_oarfish_sc_output(
        oarfish_out = x,
        annotation = annotation,
        outdir = outdir
      )
    }) |>
      setNames(samples_oarfish)

  } else {
    if (missing(quantification)) {
      stop(sprintf("No transcript count results found in %s, folder needs to contain transcript_count.csv.gz (FLAMES quantification results) or count.mtx (Oarfish quantification)", outdir))
    } else if (quantification == "FLAMES") {
      stop(sprintf("No transcript count results found in %s, folder needs to contain transcript_count.csv.gz (FLAMES quantification results)", outdir))
    } else {
      stop(sprintf("No transcript count results found in %s, folder needs to contain count.mtx (Oarfish quantification)", outdir))
    }
  }

  if (length(sce_list) == 1) {
    tryCatch({
      sce_list[[1]] <- add_gene_counts(sce_list[[1]], file.path(outdir, "gene_count.csv"))
    }, error = function(e) {
      message(sprintf("Gene counts not added to SingleCellExperiment object: %s", e$message))
    })
    return(sce_list[[1]])
  } else {
    sce_list <- sapply(names(sce_list), \(x) {
      tryCatch({
        sce_list[[x]] <- add_gene_counts(sce_list[[x]], paste(x, "gene_count.csv", sep = "_"))
      }, error = function(e) {
        message(sprintf("Gene counts not added to sample %s: %s", x, e$message))
      })
    }, simplify = FALSE)
  }

  names(sce_list) <- basename(names(sce_list))
  return(sce_list)
}

#' Add rowRanges by rownames to \code{SummarizedExperiment} object
#' Assumes rownames are transcript_ids
#' Assumes transcript_id is present in the annotation file
#' @importFrom SummarizedExperiment rowRanges rowRanges<-
#' @return a \code{SummarizedExperiment} object with rowRanges added
#' @keywords internal
addRowRanges <- function(sce, annotation, outdir) {
  if (is.null(S4Vectors::metadata(sce)$OutputFiles)) {
    S4Vectors::metadata(sce)$OutputFiles <- list()
  }

  if (file.exists(file.path(outdir, "isoform_annotated.gtf"))) {
    novel_annotation <- file.path(outdir, "isoform_annotated.gtf")
  } else if (file.exists(file.path(outdir, "isoform_annotated.gff3"))) {
    novel_annotation <- file.path(outdir, "isoform_annotated.gff3")
  } else if (!is.null(S4Vectors::metadata(sce)$OutputFiles$isoform_annotated)) {
    novel_annotation <- S4Vectors::metadata(sce)$OutputFiles$isoform_annotated
  } else {
    message(sprintf("isoform_annotated.gff3/gtf file not found in %s, this can be safely ignored if the pipeline was run with do_isoform_identification = FALSE", outdir))
    novel_annotation <- NULL
  }

  annotation_grl <- get_GRangesList(annotation)
  if (!is.null(novel_annotation)) {
    novel_grl <- get_GRangesList(novel_annotation)
    annotation_grl <- c(annotation_grl,
      novel_grl[!names(novel_grl) %in% names(annotation_grl)]
    )
    S4Vectors::metadata(sce)$OutputFiles$isoform_annotated <- novel_annotation
  }

  if (any(!rownames(sce) %in% names(annotation_grl))) {
    warning(sprintf(
      "Some transcript(s) are not recorded in the annotation file: %s",
      paste(
        head(rownames(sce)[!rownames(sce) %in% names(annotation_grl)]),
        collapse = ", "
      )
    ))
  }

  annotation_grl <- annotation_grl[names(annotation_grl) %in% rownames(sce)]

  # rowData lost when adding rowRanges
  # https://github.com/Bioconductor/SummarizedExperiment/issues/81#issuecomment-2632781880
  rowDataBackup <- SummarizedExperiment::rowData(sce)

  # SummarizedExperiment throws: invalid type/length (S4/0) in vector allocation
  # convert to RangedSummarizedExperiment
  if (!is(sce, "RangedSummarizedExperiment")) {
    # don't convert SingleCellExperiment
    stopifnot("Unexpected class when adding rowRanges" = !is(sce, "SingleCellExperiment"))
    sce <- as(sce, "RangedSummarizedExperiment")
  }
  SummarizedExperiment::rowRanges(sce)[names(annotation_grl)] <- annotation_grl
  SummarizedExperiment::rowData(sce) <- rowDataBackup
  return(sce)
}

#' Create \code{SummarizedExperiment} object from \code{FLAMES} output folder
#' @param outdir The folder containing \code{FLAMES} output files
#' @param annotation (Optional) the annotation file that was used to produce the output files
#' @param quantification (Optional)  the quantification method used to generate the
#' output files (either "FLAMES" or "Oarfish".). If not specified, the function will
#' attempt to determine the quantification method.
#' @return a \code{SummarizedExperiment} object
#' @examples
#' ppl <- example_pipeline("BulkPipeline")
#' ppl <- run_FLAMES(ppl)
#' se1 <- experiment(ppl)
#' se2 <- create_se_from_dir(ppl@outdir, ppl@annotation)
#' @export
create_se_from_dir <- function(outdir, annotation, quantification = "FLAMES") {
  if (missing(quantification)) {
    quantification <- ifelse(
      length(list.files(outdir, pattern = "^transcript_count.csv.gz$", all.files = TRUE)) > 0,
      "FLAMES",
      "Oarfish"
    )
  }
  if (quantification == "FLAMES") {
    if (missing(annotation) || !is.character(annotation)) {
      annotation_file <- NULL
    } else {
      annotation_file <- annotation
    }
    out_files <- list(
      counts = file.path(outdir, "transcript_count.csv.gz"),
      outdir = outdir,
      annotation = annotation_file,
      transcript_assembly = file.path(outdir, "transcript_assembly.fa"),
      align_bam = file.path(outdir, "align2genome.bam"),
      realign2transcript = file.path(outdir, "realign2transcript.bam"),
      tss_tes = file.path(outdir, "tss_tes.bedgraph")
    )
    se <- generate_bulk_summarized(out_files)
    if (!missing(annotation)) {
      se <- addRowRanges(se, annotation, outdir)
    }
    return(se)
  } else if (quantification == "Oarfish") {
    oarfish_samples <- list.files(outdir, pattern = "\\.quant$", full.names = TRUE, all.files = TRUE) |>
      stringr::str_remove("\\.quant$")
    se <- parse_oarfish_bulk_output(
      oarfish_outs = oarfish_samples,
      sample_names = basename(oarfish_samples)
    )
    if (!missing(annotation)) {
      se <- addRowRanges(se, annotation, outdir)
    }
    return(se)
  } else {
    stop("Unknown quantification method: ", quantification)
  }
}
