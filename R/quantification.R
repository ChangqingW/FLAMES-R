#' @importFrom reticulate import_from_path dict
#' @importFrom basilisk basiliskRun
parse_realigned_bam <- function(
    bam_in, fa_idx_f, min_sup_reads,
    min_tr_coverage, min_read_coverage, bc_file) {
  ret <- basiliskRun(
    env = flames_env,
    fun = function(bam_in, fa_idx_f, min_sup_reads,
                   min_tr_coverage, min_read_coverage, bc_file) {
      python_path <- system.file("python", package = "FLAMES")
      count <- reticulate::import_from_path("count_tr", python_path)
      ret <- count$parse_realigned_bam(
        bam_in, fa_idx_f, min_sup_reads,
        min_tr_coverage, min_read_coverage, bc_file
      )
      names(ret) <- c("bc_tr_count_dict", "bc_tr_badcov_count_dict", "tr_kept")
      ret
    },
    bam_in = bam_in, fa_idx_f = fa_idx_f, min_sup_reads = min_sup_reads,
    min_tr_coverage = min_tr_coverage, min_read_coverage = min_read_coverage,
    bc_file = bc_file
  )
  ret
}

#' @importFrom reticulate import_from_path
#' @importFrom basilisk basiliskRun
#' @importFrom future plan
wrt_tr_to_csv <- function(
    bc_tr_count_dict, transcript_dict, csv_f,
    transcript_dict_ref = NULL, has_UMI = TRUE) {
  future::plan(future::multisession)
  basiliskRun(
    env = flames_env,
    fun = function(bc_tr_count_dict, transcript_dict,
                   csv_f, transcript_dict_ref, has_UMI) {
      python_path <- system.file("python", package = "FLAMES")
      count <- reticulate::import_from_path("count_tr", python_path)
      count$wrt_tr_to_csv(bc_tr_count_dict, transcript_dict, csv_f, transcript_dict_ref, has_UMI)
    },
    bc_tr_count_dict = bc_tr_count_dict, transcript_dict = transcript_dict,
    csv_f = csv_f, transcript_dict_ref = transcript_dict_ref, has_UMI = has_UMI
  )
}


#' Gene quantification
#' @description Calculate the per gene UMI count matrix by parsing the genome alignment file.
#'
#' @details
#' After the genome alignment step (\code{do_genome_align}), the alignment file will be parsed to
#' generate the per gene UMI count matrix. For each gene in the annotation file, the number of reads
#' overlapping with the gene’s genomic coordinates will be assigned to that gene. If a read overlaps
#' multiple genes, it will be assigned to the gene with the highest number of overlapping nucleotides.
#' If exon coordinates are included in the provided annotation, the decision will first consider the
#' number of nucleotides aligned to the exons of each gene. In cases of a tie, the overlap with introns
#' will be used as a tiebreaker. If there is still a tie after considering both exons and introns,
#' a random gene will be selected from the tied candidates.
#'
#' After the read-to-gene assignment, the per gene UMI count matrix will be generated.
#' Specifically, for each gene, the reads with similar mapping coordinates of transcript
#' termination sites (TTS, i.e. the end of the the read with a polyT or polyA) will be grouped
#' together. UMIs of reads in the same group will be collapsed to generate the UMI counts for each
#' gene.
#'
#' Finally, a new fastq file with deduplicated reads by keeping the longest read in each UMI.
#'
#' @param annotation The file path to the annotation file in GFF3 format
#' @param outdir The path to directory to store all output files.
#' @param infq The input FASTQ file.
#' @param n_process The number of processes to use for parallelization.
#' @param pipeline The pipeline type as a character string, either \code{sc_single_sample} (single-cell, single-sample),
#' @param samples A vector of sample names, default to the file names of input fastq files,
#' or folder names if \code{fastqs} is a vector of folders.
#' \code{bulk} (bulk, single or multi-sample), or \code{sc_multi_sample} (single-cell, multiple samples)
#' @param random_seed The random seed for reproducibility.
#' @return The count matrix will be saved in the output folder as \code{transcript_count.csv.gz}.
#' @importFrom reticulate import_from_path dict
#' @importFrom basilisk basiliskRun
quantify_gene <- function(annotation, outdir, infq, n_process, pipeline = "sc_single_sample", samples = NULL, random_seed = 2024) {
  cat(format(Sys.time(), "%X %a %b %d %Y"), "quantify genes \n")

  if (grepl("\\.gff3?(\\.gz)?$", annotation)) {
    warning("Annotation in GFF format may cause errors. Please consider using GTF formats.\n")
  }

  genome_bam <- list.files(outdir)[grepl("_?align2genome\\.bam$", list.files(outdir))]
  cat("Found genome alignment file(s): ")
  cat(paste0("\t", paste(genome_bam, collapse = "\n\t"), "\n"))

  if (length(genome_bam) != 1 && grepl("single_sample", pipeline)) {
    stop("Incorrect number of genome alignment files found.\n")
  }

  tryCatch(
    {
      basiliskRun(
        env = flames_env, fun = function(annotation, outdir, pipeline, n_process, infq, samples, random_seed) {
          python_path <- system.file("python", package = "FLAMES")
          count <- reticulate::import_from_path("count_gene", python_path)
          count$quantification(annotation, outdir, pipeline, n_process, infq = infq, sample_names = samples, random_seed = random_seed)
        },
        annotation = annotation,
        outdir = outdir,
        pipeline = pipeline,
        n_process = n_process,
        infq = infq,
        samples = samples,
        random_seed = random_seed
      )
    },
    error = function(e) {
      # Capture the Python error using py_last_error()
      py_error <- reticulate::py_last_error()
      if (!is.null(py_error)) {
        py_error_message <- py_error$message
        # Print the actual function call
        cat(annotation, outdir, pipeline, n_process, infq, samples, random_seed)
        stop("Error when quantifying genes:\n", py_error_message)
      } else {
        stop("Error when quantifying genes:\n", e$message)
      }
    }
  )
}

#' Add gene counts to a \code{SingleCellExperiment} object
#' @description Add gene counts to a \code{SingleCellExperiment} object
#' as an \code{altExps} slot named \code{gene}.
#' @param sce A \code{SingleCellExperiment} object.
#' @param gene_count_file The file path to the gene count file. If missing,
#' the function will try to find the gene count file in the output directory.
#' @return A \code{SingleCellExperiment} object with gene counts added.
#' @importFrom SingleCellExperiment SingleCellExperiment altExps
#' @importFrom S4Vectors metadata
#' @examples
#' # Set up a mock SingleCellExperiment object
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(counts = matrix(0, nrow = 10, ncol = 10))
#' )
#' colnames(sce) <- paste0("cell", 1:10)
#' # Set up a mock gene count file
#' gene_count_file <- tempfile()
#' gene_mtx <- matrix(1:10, nrow = 2, ncol = 5)
#' colnames(gene_mtx) <- paste0("cell", 1:5)
#' rownames(gene_mtx) <- c("gene1", "gene2")
#' write.csv(gene_mtx, gene_count_file)
#' # Add gene counts to the SingleCellExperiment object
#' sce <- add_gene_counts(sce, gene_count_file)
#' # verify the gene counts are added
#' SingleCellExperiment::altExps(sce)$gene
#' @export
add_gene_counts <- function(sce, gene_count_file) {
  if (missing(gene_count_file)) {
    tryCatch(
      {
        gene_count_file <- file.path(
          S4Vectors::metadata(sce)$inputs$outdir,
          "gene_count.csv"
        )
      },
      error = function(e) {
        stop("Error when finding gene count file:\n", e$message)
      }
    )
    if (!file.exists(gene_count_file)) {
      stop("Gene count file not found.")
    }
  }
  mtx <- read.csv(gene_count_file, row.names = 1)
  mtx <- mtx[, colnames(mtx) %in% colnames(sce)] |>
    as.matrix()
  gene_mtx <- matrix(0, nrow = nrow(mtx), ncol = ncol(sce))
  gene_mtx[, match(colnames(mtx), colnames(sce))] <- mtx

  gene <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = gene_mtx)
  )
  colnames(gene) <- colnames(sce)
  rownames(gene) <- rownames(mtx)

  SingleCellExperiment::altExps(sce)$gene <- gene
  return(sce)
}

#' FLAMES Transcript quantification
#' @description Calculate the transcript count matrix by parsing the re-alignment file.
#' @param annotation The file path to the annotation file in GFF3 format
#' @param outdir The path to directory to store all output files.
#' @param config Parsed FLAMES configurations.
#' @param pipeline The pipeline type as a character string, either \code{sc_single_sample} (single-cell, single-sample),
#' @param samples A vector of sample names, required for \code{sc_multi_sample} pipeline.
#' \code{bulk} (bulk, single or multi-sample), or \code{sc_multi_sample} (single-cell, multiple samples)
#' @return A \code{SingleCellExperiment} object for single-cell pipeline, a list of \code{SingleCellExperiment} objects for multi-sample pipeline, or a \code{SummarizedExperiment} object for bulk pipeline.
#' @importFrom basilisk basiliskRun
#' @importFrom reticulate import_from_path dict
#' @keywords internal
quantify_transcript_flames <- function(annotation, outdir, config, pipeline = "sc_single_sample", samples) {
  cat(format(Sys.time(), "%X %a %b %d %Y"), "quantify transcripts \n")

  if (grepl("\\.gff3?(\\.gz)?$", annotation)) {
    warning("Annotation in GFF format may cause errors. Please consider using GTF formats.\n")
  }

  realign_bam <- list.files(outdir)[grepl("_?realign2transcript\\.bam$", list.files(outdir))]
  cat("Found realignment file(s): ")
  cat(paste0("\t", paste(realign_bam, collapse = "\n\t"), "\n"))

  if (length(realign_bam) != 1 && grepl("single_sample", pipeline)) {
    stop("Incorrect number of realignment files found.\n")
  }

  known_transcripts <- annotation |>
    rtracklayer::import(feature.type = "transcript") |>
    S4Vectors::mcols() |>
    (\(x) x$transcript_id)()

  basiliskRun(
    env = flames_env, fun = function(config_dict, annotation, known_transcripts, outdir, pipeline) {
      python_path <- system.file("python", package = "FLAMES")
      count <- reticulate::import_from_path("count_tr", python_path)
      count$quantification(config_dict, annotation, known_transcripts, outdir, pipeline)
    },
    config_dict = reticulate::dict(config),
    annotation = annotation,
    known_transcripts = known_transcripts,
    outdir = outdir,
    pipeline = pipeline
  )

  if (pipeline == "sc_single_sample") {
    out_files <- list(
      "annotation" = annotation,
      # "genome_fa" = genome_fa,
      "counts" = file.path(outdir, "transcript_count.csv.gz"),
      "isoform_annotated" = file.path(outdir, "isoform_annotated.filtered.gff3"),
      "transcript_assembly" = file.path(outdir, "transcript_assembly.fa"),
      # "align_bam" = genome_bam,
      "realign2transcript" = file.path(outdir, "realign2transcript.bam"),
      # "tss_tes" = file.path(outdir, "tss_tes.bedgraph"),
      "fsm_annotation" = file.path(outdir, "isoform_FSM_annotation.csv"),
      "outdir" = outdir
    )
    sce <- generate_sc_singlecell(out_files) |>
      addRowRanges(annotation, outdir)
    return(sce)
  } else if (pipeline == "sc_multi_sample") {
    sce_list <- as.list(seq_along(samples))
    names(sce_list) <- samples
    for (i in seq_along(samples)) {
      out_files <- list(
        "annotation" = annotation,
        "counts" = file.path(outdir, paste0(samples[i], "_transcript_count.csv.gz")),
        "isoform_annotated" = file.path(outdir, ifelse(config$pipeline_parameters$bambu_isoform_identification, "isoform_annotated.gtf", "isoform_annotated.gff3")),
        "transcript_assembly" = file.path(outdir, "transcript_assembly.fa"),
        "realign2transcript" = file.path(outdir, paste0(samples[i], "_realign2transcript.bam")),
        "outdir" = outdir,
        "fsm_annotation" = file.path(outdir, "isoform_FSM_annotation.csv")
      )
      sce_list[[i]] <- generate_sc_singlecell(out_files) |>
        addRowRanges(annotation, outdir)
    }
    return(sce_list)
  } else if (pipeline == "bulk") {
    out_files <- list(
      "annotation" = annotation,
      "counts" = file.path(outdir, "transcript_count.csv.gz"),
      "isoform_annotated" = file.path(outdir, "isoform_annotated.filtered.gff3"),
      "transcript_assembly" = file.path(outdir, "transcript_assembly.fa"),
      "realign2transcript" = file.path(outdir, list.files(outdir))[grepl("realign2transcript\\.bam$", list.files(outdir))],
      "outdir" = outdir,
      "fsm_annotation" = file.path(outdir, "isoform_FSM_annotation.csv")
    )
    se <- generate_bulk_summarized(out_files) |>
      addRowRanges(annotation, outdir)
    return(se)
  }
}

#' @importFrom basilisk obtainEnvironmentPath basiliskRun
run_oarfish <- function(realign_bam, outdir, threads = 1, sample, oarfish_bin, single_cell = TRUE) {
  if (missing(oarfish_bin)) {
    oarfish_bin <- find_bin("oarfish")
    stopifnot(!is.na(oarfish_bin))
  }

  if (missing(sample)) {
    sample <- "oarfish"
  }

  oarfish_status <- base::system2(
    command = oarfish_bin,
    args = c(
      switch(single_cell,
        "--single-cell"
      ),
      "--alignments", file.path(outdir, realign_bam),
      "-j", threads, "--output", file.path(outdir, sample)
    ),
  )
  if (oarfish_status != 0) {
    stop(paste0("error running oarfish:\n", oarfish_status))
  }

  return(file.path(outdir, sample))
}

#' @importFrom MatrixGenerics rowSums
#' @importFrom Matrix readMM
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment rowRanges rowRanges<-
parse_oarfish_sc_output <- function(oarfish_out, annotation, outdir) {
  mtx <- t(Matrix::readMM(paste0(oarfish_out, ".count.mtx")))
  rownames(mtx) <- read.delim(paste0(oarfish_out, ".features.txt"), header = FALSE)$V1
  colnames(mtx) <- read.delim(paste0(oarfish_out, ".barcodes.txt"), header = FALSE)$V1
  # mtx <- mtx[MatrixGenerics::rowSums(mtx) > 0, ]
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mtx))

  annotation_grl <- get_GRangesList(annotation)
  if (file.exists(file.path(outdir, "isoform_annotated.gff3"))) {
    novel_annotation <- file.path(outdir, "isoform_annotated.gff3")
  } else if (file.exists(file.path(outdir, "isoform_annotated.gtf"))) {
    novel_annotation <- file.path(outdir, "isoform_annotated.gtf")
  } else {
    novel_annotation <- NULL
    message("No novel annotation found.")
  }
  if (!is.null(novel_annotation)) {
    novel_grl <- get_GRangesList(novel_annotation)
    annotation_grl <- c(
      annotation_grl,
      novel_grl[!names(novel_grl) %in% names(annotation_grl)]
    )
  }

  annotation_grl <- annotation_grl[names(annotation_grl) %in% rownames(sce)]
  SummarizedExperiment::rowRanges(sce)[names(annotation_grl)] <- annotation_grl
  rowData(sce)$transcript_id <- rownames(sce)

  # add gene ID
  SummarizedExperiment::rowData(sce)$gene_id <- tryCatch(
    GenomicRanges::mcols(annotation_grl)[, "gene_id"],
    error = function(e) {
      sprintf("Error when adding gene ID: %s, trying again with txdbmaker...", e$message)
      gene_id_tb <- c(annotation, novel_annotation) |>
        lapply(fake_stranded_gff) |>
        lapply(txdbmaker::makeTxDbFromGFF) |>
        lapply(\(x) GenomicFeatures::transcripts(x, columns = c("gene_id", "tx_name"))) |>
        lapply(\(x) GenomicRanges::mcols(x)) |>
        do.call(rbind, args = _)
      gene_id_tb[match(rownames(sce), gene_id_tb$tx_name), "gene_id"] |>
        as.character()
    }
  )

  return(sce)
}

#' @importFrom SummarizedExperiment SummarizedExperiment
parse_oarfish_bulk_output <- function(oarfish_outs, sample_names) {
  mtx_list <- lapply(oarfish_outs, function(oarfish_out) {
    read.delim(paste0(oarfish_out, ".quant"), header = TRUE, row.names = 1)[, "num_reads", drop = FALSE]
  })
  mtx <- do.call(cbind, mtx_list)
  colnames(mtx) <- sample_names
  SummarizedExperiment::SummarizedExperiment(assays = list(counts = mtx))
}

quantify_transcript_oarfish <- function(
    annotation, outdir, config,
    pipeline = "sc_single_sample", samples) {
  realign_bam <- list.files(outdir)[grepl("_?realign2transcript\\.bam$", list.files(outdir))]
  if (pipeline == "sc_single_sample") {
    oarfish_out <- run_oarfish(realign_bam, outdir, threads = config$pipeline_parameters$threads)
    return(parse_oarfish_sc_output(oarfish_out, annotation, outdir))
  } else if (pipeline == "sc_multi_sample") {
    sce_list <- as.list(1:length(samples))
    names(sce_list) <- samples
    for (i in 1:length(samples)) {
      realign_bam <- list.files(outdir)[grepl(paste0(samples[i], "_realign2transcript\\.bam$"), list.files(outdir))]
      oarfish_out <- run_oarfish(realign_bam, outdir, threads = config$pipeline_parameters$threads, sample = samples[i])
      sce_list[[i]] <- parse_oarfish_sc_output(oarfish_out, annotation, outdir)
    }
    return(sce_list)
  } else if (pipeline == "bulk") {
    oarfish_out <- rep(NA, length(samples))
    for (i in 1:length(samples)) {
      realign_bam <- list.files(outdir)[grepl(paste0(samples[i], "_realign2transcript\\.bam$"), list.files(outdir))]
      oarfish_out[i] <- run_oarfish(realign_bam, outdir, threads = config$pipeline_parameters$threads, sample = samples[i], single_cell = FALSE)
    }
    return(parse_oarfish_bulk_output(oarfish_out, samples))
  } else {
    stop(paste0("Unknown pipeline: ", pipeline))
  }
}

#' Transcript quantification
#' @description Calculate the transcript count matrix by parsing the re-alignment file.
#' @param annotation The file path to the annotation file in GFF3 format
#' @param outdir The path to directory to store all output files.
#' @param config Parsed FLAMES configurations.
#' @param pipeline The pipeline type as a character string, either \code{sc_single_sample} (single-cell, single-sample),
#' @param ... Supply sample names as character vector (e.g. \code{samples = c("name1", "name2", ...)}) for muti-sample or bulk pipeline.
#' \code{bulk} (bulk, single or multi-sample), or \code{sc_multi_sample} (single-cell, multiple samples)
#' @return A \code{SingleCellExperiment} object for single-cell pipeline, a list of \code{SingleCellExperiment} objects for multi-sample pipeline, or a \code{SummarizedExperiment} object for bulk pipeline.
#' @keywords internal
quantify_transcript <- function(annotation, outdir, config, pipeline = "sc_single_sample", ...) {
  if (config$pipeline_parameters$oarfish_quantification) {
    res <- quantify_transcript_oarfish(annotation, outdir, config, pipeline, ...)
  } else {
    res <- quantify_transcript_flames(annotation, outdir, config, pipeline, ...)
  }
  if (is.list(res)) {
    # sce <- do.call(BiocGenerics::cbind, res)
    # colData(sce)$sample <- unlist(mapply(function(x, n) rep(x, n), names(res), sapply(res, ncol)))
    return(res)
  } else {
    return(res)
  }
}

# example for Rsamtools
#' @importFrom Rsamtools BamFile  scanBam isIncomplete
# quantify_tmp <- function(bamFileName) {
#    bf <- Rsamtools::BamFile(bamFileName, yieldSize=100)
#    while (Rsamtools::isIncomplete(bf)) {
#        print(Rsamtools::scanBam(bf)[[1]][[1]])
#        Sys.sleep(3)
#    }
# }
