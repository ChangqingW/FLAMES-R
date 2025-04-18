#' Pipeline for Single Cell Data
#'
#' @md
#'
#' @description
#' Semi-supervised isoform detection and annotation for long read data.
#' This variant is for single cell data. By default, this pipeline demultiplexes input
#' fastq data (\code{match_cell_barcode = TRUE}). Specific parameters relating to
#' analysis can be changed either through function arguments, or through a
#' configuration JSON file.
#'
#' @details
#' By default FLAMES use minimap2 for read alignment. After the genome alignment step (\code{do_genome_align}), FLAMES summarizes the alignment for each read by grouping reads
#' with similar splice junctions to get a raw isoform annotation (\code{do_isoform_id}). The raw isoform
#' annotation is compared against the reference annotation to correct potential splice site
#' and transcript start/end errors. Transcripts that have similar splice junctions
#' and transcript start/end to the reference transcript are merged with the
#' reference. This process will also collapse isoforms that are likely to be truncated
#' transcripts. If \code{isoform_id_bambu} is set to \code{TRUE}, \code{bambu::bambu} will be used to generate the updated annotations.
#' Next is the read realignment step (\code{do_read_realign}), where the sequence of each transcript from the update annotation is extracted, and
#' the reads are realigned to this updated \code{transcript_assembly.fa} by minimap2. The
#' transcripts with only a few full-length aligned reads are discarded.
#' The reads are assigned to transcripts based on both alignment score, fractions of
#' reads aligned and transcript coverage. Reads that cannot be uniquely assigned to
#' transcripts or have low transcript coverage are discarded. The UMI transcript
#' count matrix is generated by collapsing the reads with the same UMI in a similar
#' way to what is done for short-read scRNA-seq data, but allowing for an edit distance
#' of up to 2 by default. Most of the parameters, such as the minimal distance to splice site and minimal percentage of transcript coverage
#' can be modified by the JSON configuration file (\code{config_file}).
#'
#' @param annotation The file path to the annotation file in GFF3 format
#' @param fastq The file path to input fastq file
#' @param genome_bam Optional file path to a bam file to use instead of fastq file (skips initial alignment step)
#' @param outdir The path to directory to store all output files.
#' @param genome_fa The file path to genome fasta file.
#' @param minimap2 Path to minimap2, if it is not in PATH. Only required if either or both of
#' \code{do_genome_align} and \code{do_read_realign} are \code{TRUE}.
#' @param k8 Path to the k8 Javascript shell binary. Only required if \code{do_genome_align} is \code{TRUE}.
#' @param config_file File path to the JSON configuration file. If specified, \code{config_file} overrides
#' all configuration parameters
#' @param barcodes_file The file path to the reference csv used for demultiplexing in flexiplex. If not specified, the
#'                           demultiplexing will be performed using BLAZE. Default is \code{NULL}.
#' @param expect_cell_number Expected number of cells for identifying the barcode list in BLAZE.
#'                           This could be just a rough estimate. E.g., the targeted number of cells.
#'                           Required if the \code{do_barcode_demultiplex} are \code{TRUE} in the the JSON configuration file
#'                           and \code{barcodes_file} is not specified. Default is \code{NULL}.
#'
#' @return if \code{do_transcript_quantification} set to true, \code{sc_long_pipeline} returns a \code{SingleCellExperiment} object, containing a count
#' matrix as an assay, gene annotations under metadata, as well as a list of the other
#' output files generated by the pipeline. The pipeline also outputs a number of output
#' files into the given \code{outdir} directory. These output files generated by the pipeline are:
#' \describe{
#'  \item{transcript_count.csv.gz}{ - a transcript count matrix (also contained in the SingleCellExperiment)}
#'  \item{isoform_annotated.filtered.gff3}{ - isoforms in gff3 format (also contained in the SingleCellExperiment)}
#'  \item{transcript_assembly.fa}{ - transcript sequence from the isoforms}
#'  \item{align2genome.bam}{ - sorted BAM file with reads aligned to genome}
#'  \item{realign2transcript.bam}{ - sorted realigned BAM file using the transcript_assembly.fa as reference}
#'  \item{tss_tes.bedgraph}{ - TSS TES enrichment for all reads (for QC)}
#' }
#' if \code{do_transcript_quantification} set to false, nothing will be returned
#'
#' @details The default parameters can be changed either through the function
#' arguments are through the configuration JSON file \code{config_file}. the \code{pipeline_parameters}
#' section specifies which steps are to be executed in the pipeline - by default, all
#' steps are executed. The \code{isoform_parameters} section affects isoform detection - key
#' parameters include:
#' \describe{
#'  \item{\code{Min_sup_cnt}}{ which causes transcripts with less reads aligned than
#' it's value to be discarded}
#'  \item{\code{MAX_TS_DIST}}{ which merges transcripts with the same intron
#' chain and TSS/TES distace less than \code{MAX_TS_DIST}}
#'  \item{\code{strand_specific}}{ which specifies if reads are in the same strand as the mRNA (1),
#' or the reverse complemented (-1) or not strand specific (0), which results in
#' strand information being based on reference annotation.}
#' }
#'
#' @seealso
#' [bulk_long_pipeline()] for bulk long data,
#' [SingleCellExperiment()] for how data is outputted
#'
#' @importFrom dplyr group_by summarise_at slice_max filter
#' @importFrom magrittr "%>%"
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDimNames logcounts
#' @importFrom SummarizedExperiment rowData colData rowData<- colData<- rowRanges rowRanges<-
#' @importFrom BiocGenerics cbind colnames rownames start end
#' @importFrom utils read.csv read.table
#' @importFrom GenomeInfoDb seqlengths
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
#' if (!any(is.na(find_bin(c("minimap2", "k8"))))) {
#'   sce <- FLAMES::sc_long_pipeline(
#'     genome_fa = genome_fa,
#'     fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
#'     annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
#'     outdir = outdir,
#'     barcodes_file = bc_allow
#'   )
#' }
#' @export
sc_long_pipeline <- function(
    annotation, fastq, genome_bam = NULL, outdir, genome_fa,
    minimap2 = NULL, k8 = NULL, barcodes_file = NULL, expect_cell_number = NULL, config_file = NULL) {

  # Initialize the pipeline object
  pipeline <- new("FLAMES.SingleCellPipeline")
  config <- check_arguments(annotation, fastq, genome_bam, outdir, genome_fa, config_file)$config
  steps <- c("demultiplexing", "alignment", "quantification", "isoform_identification", "realignment")
  pipeline <- initializePipeline(pipeline, config, steps, outdir)

  # Assign specific slots
  pipeline@fastq <- fastq
  pipeline@genome_bam <- genome_bam
  pipeline@barcodes_file <- barcodes_file
  pipeline@expect_cell_number <- expect_cell_number

  # Example: Accessing pipeline slots
  cat("Pipeline initialized with steps:", paste(pipeline@steps, collapse = ", "), "\n")

  # Continue with the existing logic, refactored to use the pipeline object
  metadata <- list(
    "inputs" = list(
      "config_file" = config_file,
      "annotation" = annotation,
      "fastq" = fastq,
      "outdir" = outdir,
      "genome_fa" = genome_fa,
      "barcodes_file" = barcodes_file,
      "expect_cell_number" = expect_cell_number,
      "minimap2" = minimap2,
      "k8" = k8,
      "genome_bam" = genome_bam
    ),
    "results" = list()
  )
  metadata$inputs <- metadata$inputs[!sapply(metadata$inputs, is.null)]

  # Demultiplexing step
  if (config$pipeline_parameters$do_barcode_demultiplex) {
    if (is.null(pipeline@barcodes_file)) {
      cat("Running BLAZE to generate barcode list from long reads...\n")
      if (is.null(pipeline@expect_cell_number)) {
        stop("'expect_cell_number' is required to run BLAZE for barcode identification. Please specify it.")
      }
      blaze(pipeline@expect_cell_number, pipeline@fastq,
        "output-prefix" = paste0(outdir, "/"),
        "output-fastq" = "matched_reads.fastq",
        "threads" = config$pipeline_parameters$threads,
        "max-edit-distance" = config$barcode_parameters$max_bc_editdistance,
        "overwrite" = TRUE
      )
      pipeline@fastq <- file.path(outdir, "matched_reads.fastq")
    } else {
      cat("Demultiplexing using flexiplex...\n")
      pipeline@fastq <- file.path(outdir, "matched_reads.fastq")
      bc_stat <- file.path(outdir, "matched_barcode_stat")
      metadata$results <- c(metadata$results,
        list("find_barcode" = find_barcode(
          fastq = pipeline@fastq,
          barcodes_file = pipeline@barcodes_file,
          stats_out = bc_stat,
          reads_out = pipeline@fastq,
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
        )
        )
      )
    }
  }

  # Align reads to genome
  if (config$pipeline_parameters$do_genome_alignment) {
    cat("#### Aligning reads to genome using minimap2\n")
    metadata$results <- c(metadata$results,
      list("minimap2_align" = minimap2_align(config, genome_fa, pipeline@fastq,
        annotation, outdir, minimap2, k8, prefix = NULL,
        threads = config$pipeline_parameters$threads
      )
      )
    )
  }

  # Gene quantification and UMI deduplication
  if (config$pipeline_parameters$do_gene_quantification) {
    cat("#### Start gene quantification and UMI deduplication\n")
    quantify_gene(annotation, outdir, pipeline@fastq, config$pipeline_parameters$threads,
      pipeline = "sc_single_sample", random_seed = config$pipeline_parameters$seed
    )
  }

  # Isoform identification
  if (config$pipeline_parameters$do_isoform_identification) {
    find_isoform(annotation, genome_fa, pipeline@genome_bam, outdir, config)
  }

  # Realign to transcript
  if (config$pipeline_parameters$do_read_realignment) {
    cat("#### Realigning deduplicated reads to transcript using minimap2\n")
    infq_realign <- if (config$pipeline_parameters$do_gene_quantification) {
      file.path(outdir, "matched_reads_dedup.fastq")
    } else {
      pipeline@fastq
    }
    metadata$results <- c(metadata$results,
      list("minimap2_realign" = minimap2_realign(config, infq_realign, outdir, minimap2,
        prefix = NULL, threads = config$pipeline_parameters$threads)
      )
    )
  }

  # Transcript quantification
  if (config$pipeline_parameters$do_transcript_quantification) {
    sce <- quantify_transcript(
      annotation = annotation,
      outdir = outdir,
      pipeline = "sc_single_sample",
      config = config
    )
    sce@metadata <- c(sce@metadata, metadata)
    return(sce)
  }

  return(metadata)
}

#' @importFrom utils read.csv
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
#' sce <- FLAMES::sc_long_pipeline(
#'   genome_fa = genome_fa,
#'   fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
#'   annotation = annotation,
#'   outdir = outdir,
#'   barcodes_file = bc_allow,
#'   config_file = create_config(outdir, oarfish_quantification = FALSE)
#' )
#' sce_2 <- create_sce_from_dir(outdir, annotation)
create_sce_from_dir <- function(outdir, annotation, quantification = "FLAMES") {

  samples <- list.files(outdir, pattern = "_?transcript_count.csv.gz$", full.names = TRUE)
  samples_oarfish <- list.files(outdir, pattern = "\\.count\\.mtx$", full.names = TRUE) |>
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
        addRowRanges(annotation, outdir) |>
        setNames(stringr::str_remove(x, "_?transcript_count.csv.gz$"))
    })

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
  SummarizedExperiment::rowRanges(sce)[names(annotation_grl)] <- annotation_grl
  SummarizedExperiment::rowData(sce) <- rowDataBackup
  return(sce)
}

#' Create \code{SummarizedExperiment} object from \code{FLAMES} output folder
#' @param outdir The folder containing \code{FLAMES} output files
#' @param annotation (Optional) the annotation file that was used to produce the output files
#' @return a \code{SummarizedExperiment} object
#' @example inst/examples/pipeline_example.R
#' @export
create_se_from_dir <- function(outdir, annotation) {
  out_files <- list(
    counts = file.path(outdir, "transcript_count.csv.gz"),
    outdir = outdir,
    annotation = annotation,
    transcript_assembly = file.path(outdir, "transcript_assembly.fa"),
    align_bam = file.path(outdir, "align2genome.bam"),
    realign2transcript = file.path(outdir, "realign2transcript.bam"),
    tss_tes = file.path(outdir, "tss_tes.bedgraph")
  )
  se <- generate_bulk_summarized(out_files) |>
    addRowRanges(annotation, outdir)
  return(se)
}
