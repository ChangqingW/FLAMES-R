# Classes for segment and group
setClass("FlexiplexSegment",
         slots = list(
           type = "character",
           pattern = "character",
           name = "character",
           bc_list_name = "character",
           group = "character",
           buffer_size = "numeric",
           max_edit_distance = "numeric"
         ),
         prototype = list(
           type = "FIXED",
           bc_list_name = NA_character_,
           group = NA_character_,
           buffer_size = 2,
           max_edit_distance = 2
         ))
barcode_segment <- function(type = "FIXED", pattern, name,
                    bc_list = NA_character_, group = NA_character_,
                    buffer_size = 2, max_edit_distance = 2) {
  if (!(type %in% c("FIXED", "MATCHED", "MATCHED_SPLIT", "RANDOM"))) {
    stop("type must be one of 'FIXED', 'MATCHED', 'MATCHED_SPLIT', or 'RANDOM'")
  }

  if (!is.character(pattern) || length(pattern) != 1 || is.na(pattern) || nchar(pattern) == 0) {
    stop("pattern must be a non-empty character scalar")
  }

  if (missing(name) && type == "FIXED") {
      name <- "FIXED_SEGMENT"
  } else if (!is.character(name) || length(name) != 1 || is.na(name) || nchar(name) == 0) {
    stop("name must be a non-empty character scalar")
  }

  if (type %in% c("MATCHED", "MATCHED_SPLIT")) {
    if (!is.numeric(buffer_size) || length(buffer_size) != 1 || is.na(buffer_size)) {
      stop("buffer_size must be a numeric scalar")
    }
    if (buffer_size < 0) {
      stop("buffer_size must be non-negative")
    }
    if (!is.numeric(max_edit_distance) || length(max_edit_distance) != 1 || is.na(max_edit_distance)) {
      stop("max_edit_distance must be a numeric scalar")
    }
    if (max_edit_distance < 0) {
      stop("max_edit_distance must be non-negative")
    }
  }

  if (type == "MATCHED") {
    if (!is.character(bc_list) || length(bc_list) != 1 || is.na(bc_list) || nchar(bc_list) == 0) {
      stop("bc_list must be provided as a non-empty character scalar for type 'MATCHED'")
    }
  }
  if (type == "MATCHED_SPLIT") {
    if (!is.character(group) || length(group) != 1 || is.na(group) || nchar(group) == 0) {
      stop("group must be provided as a non-empty character scalar for type 'MATCHED_SPLIT'")
    }
  }

  new("FlexiplexSegment",
      type = type,
      pattern = pattern,
      name = name,
      bc_list_name = bc_list,
      group = group,
      buffer_size = buffer_size,
      max_edit_distance = max_edit_distance 
  )
}

setClass("FlexiplexGroup",
         slots = list(
           name = "character",
           bc_list_name = "character",
           max_edit_distance = "numeric"
         ))
barcode_group <- function(name, bc_list_name, max_edit_distance = 2) {
  if (!is.character(name) || length(name) != 1 || is.na(name) || nchar(name) == 0) {
    stop("name must be a non-empty character scalar")
  }
  if (!is.character(bc_list_name) || length(bc_list_name) != 1 || is.na(bc_list_name) || nchar(bc_list_name) == 0) {
    stop("bc_list_name must be a non-empty character scalar")
  }
  if (!is.numeric(max_edit_distance) || length(max_edit_distance) != 1 || is.na(max_edit_distance)) {
    stop("max_edit_distance must be a numeric scalar")
  }
  if (max_edit_distance < 0) {
    stop("max_edit_distance must be non-negative")
  }
  new("FlexiplexGroup",
      name = name,
      bc_list_name = bc_list_name,
      max_edit_distance = max_edit_distance)
}

#' Convert legacy pattern to Flexiplex segments
#' @param pattern named character vector defining the barcode pattern
#' @param barcodes_file path to file containing barcode allow-list
#' @param max_bc_editdistance max edit distances for the barcode sequence
#' @param buffer_size buffer size for barcode matching
#' @return a list of FlexiplexSegment objects
.legacy_pattern_to_flexiplex_segments <- function(
    pattern,
    barcodes_file,
    max_bc_editdistance,
    buffer_size
) {
  message("Converting legacy `pattern` argument to `segments`...")
  if (!is.character(pattern) || length(pattern) < 1) {
    stop("pattern must be a named character vector")
  }
  if (is.null(names(pattern)) || any(names(pattern) == "")) {
    stop("pattern must be a named character vector")
  }
  if (!"BC" %in% names(pattern)) {
    stop("pattern must contain a 'BC' entry")
  }
  if (length(barcodes_file) != 1 || !is.character(barcodes_file) || is.na(barcodes_file) || nchar(barcodes_file) == 0) {
    stop("barcodes_file must be a single non-empty string")
  }

  segments <- list()
  for (i in seq_along(pattern)) {
    part_name <- names(pattern)[i]
    part_pattern <- pattern[[i]]
    if (part_name == "BC") {
      segment <- barcode_segment(
        type = "MATCHED",
        pattern = part_pattern,
        name = "CB",
        bc_list = barcodes_file,
        buffer_size = buffer_size,
        max_edit_distance = max_bc_editdistance
      )
    } else if (part_name == "UMI") {
      segment <- barcode_segment(
        type = "RANDOM",
        pattern = part_pattern,
        name = "UB"
      )
    } else {
      segment <- barcode_segment(
        type = "FIXED",
        pattern = part_pattern,
        name = part_name
      )
    }
    segments[[length(segments) + 1]] <- segment
  }
  return(segments)
}

#' Resolve bc_list_name keys/indices to file paths
#' @param x list of FlexiplexSegment objects or list of FlexiplexGroup objects
#' @param type character, either "segment" or "barcode_group"
#' @param barcodes_files named or unnamed (only allowd when length is 1) character vector of barcode file paths
.resolve_bc_lists <- function(x, barcodes_files) {
  if (!is.null(barcodes_files)) {
    if (!is.character(barcodes_files)) {
      stop("`barcodes_files` must be a character vector of file paths")
    }
  }
  if (length(barcodes_files) > 1 && is.null(names(barcodes_files))) {
    stop("`barcodes_files` must be a named character vector when length > 1")
  }

  for (i in seq_along(x)) {
    s <- x[[i]]
    if (class(s) == "FlexiplexSegment") {
      class <- "segment"
    } else if (class(s) == "FlexiplexGroup") {
      class <- "barcode_group"
    } else {
      stop("`x` must be a list of FlexiplexSegment or FlexiplexGroup objects")
    }

    if (class == "segment" && s@type != "MATCHED") {
      next
    }

    # Single unnamed entry: use it directly
    if (length(barcodes_files) == 1 && is.null(names(barcodes_files))) {
      s@bc_list_name <- unname(barcodes_files[[1]])
      x[[i]] <- s
      next
    }

    # Get the key to resolve
    key <- s@bc_list_name
    if (length(key) != 1 || is.na(key) || nchar(key) == 0) {
      stop(sprintf(
        "%s[[%d]]: bc_list_name is empty, must provide a valid name from `barcodes_files` for barcode groups and MATCHED segments."
      ))
    }

    # Resolve via named lookup
    if (!is.null(names(barcodes_files)) && key %in% names(barcodes_files)) {
      s@bc_list_name <- unname(barcodes_files[[key]])
      x[[i]] <- s
      next
    }

    # If key directly points to an existing file, accept it.
    if (file.exists(key)) {
      message(
        sprintf(
          "%s[[%d]]: bc_list_name='%s' does not match any name in `barcodes_files`, using it as a file path directly.",
          class, i, key
        )
      )
      next
    }

    stop(
      sprintf(
        "%s[[%d]]: bc_list_name='%s' does not match any name in `barcodes_files`.", 
        class, i, key
      )
    )
  }

 x 
}

#' Match Cell Barcodes
#'
#' @description demultiplex reads with flexiplex
#' @details
#' This function demultiplexes reads by searching for flanking sequences (adaptors)
#' around the barcode sequence, and then matching against allowed barcodes.
#'
#' @importFrom utils file_test
#' @importFrom readr read_delim col_character col_integer col_logical
#' @importFrom dplyr bind_rows
#' @param fastq A path to a FASTQ file or a directory containing FASTQ files.
#' @param barcodes_file path to file containing barcode allow-list, with one barcode in each line
#' @param max_bc_editdistance max edit distances for the barcode sequence
#' @param max_flank_editdistance max edit distances for the flanking sequences (primer and polyT)
#' @param reads_out path to output FASTQ file
#' @param stats_out path of output stats file
#' @param threads number of threads to be used
#' @param cutadapt_minimum_length minimum read length after TSO trimming (cutadapt's --minimum-length)
#' @param full_length_only boolean, when TSO sequence is provided, whether reads without TSO
#' are to be discarded
#' @param pattern named character vector defining the barcode pattern
#' @param strand strand of the barcode pattern, either '+' or '-' (read will be reverse complemented
#' after barcode matching if '-')
#' @param TSO_seq TSO sequence to be trimmed
#' @param TSO_prime either 3 (when \code{TSO_seq} is on 3' the end) or 5 (on 5' end)
#' @return a list containing: \code{reads_tb} (tibble of read demultiplexed information) and
#'  \code{input}, \code{output}, \code{read1_with_adapter} from cutadapt report
#'  (if TSO trimming is performed)
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' bc_allow <- file.path(outdir, "bc_allow.tsv")
#' R.utils::gunzip(
#'   filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
#'   destname = bc_allow, remove = FALSE
#' )
#' find_barcode(
#'   fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
#'   stats_out = file.path(outdir, "bc_stat.tsv.gz"),
#'   reads_out = file.path(outdir, "demultiplexed.fastq.gz"),
#'   barcodes_file = bc_allow, 
#'   TSO_seq = "AAGCAGTGGTATCAACGCAGAGTACATGGG", TSO_prime = 5,
#'   strand = '-', cutadapt_minimum_length = 10, full_length_only = TRUE
#' )
#' @md
#' @export
find_barcode <- function(
    fastq,
    segments, barcode_groups,
    barcodes_files,
    max_flank_editdistance = 8,
    reads_out, stats_out, threads = 1,
    TSO_seq = "", TSO_prime = 3, strand = '+', cutadapt_minimum_length = 1, full_length_only = FALSE,
    sample_names = NULL,
    # legacy style
    pattern = c(
      primer = "CTACACGACGCTCTTCCGATCT",
      BC = paste0(rep("N", 16), collapse = ""),
      UMI = paste0(rep("N", 12), collapse = ""),
      polyT = paste0(rep("T", 9), collapse = "")
    ), max_bc_editdistance = 2
) {

  if (strand == '+') {
    reverseCompliment <- FALSE
  } else if (strand == '-') {
    reverseCompliment <- TRUE
  } else {
    stop("strand must be either '+' or '-'")
  }

  # Backward compatibility for legacy `pattern` arguments
  if (missing(segments) || is.null(segments) || length(segments) == 0) {
    segments <- .legacy_pattern_to_flexiplex_segments(
      pattern = pattern,
      barcodes_file = barcodes_files,
      max_bc_editdistance = max_bc_editdistance,
      buffer_size = 5
    )
    barcode_groups <- list()
  } else {
    if (!missing(pattern) && !is.null(pattern) && length(pattern) > 0) {
      warning("Both `segments` and legacy `pattern` arguments are provided, ignoring `pattern`.")
    }
    if (missing(barcode_groups)) {
      barcode_groups <- list()
    }

    if (is.data.frame(segments)) {
      segments <- apply(segments, 1, function(row) {
        barcode_segment(
          type = row[["type"]],
          pattern = row[["pattern"]],
          name = row[["name"]],
          bc_list = ifelse("bc_list_name" %in% names(row), row[["bc_list_name"]], NA_character_),
          group = ifelse("group" %in% names(row), row[["group"]], NA_character_),
          buffer_size = ifelse("buffer_size" %in% names(row), as.numeric(row[["buffer_size"]]), 2),
          max_edit_distance = ifelse("max_edit_distance" %in% names(row), as.numeric(row[["max_edit_distance"]]), 2)
        )
      })
    }
    if (is.data.frame(barcode_groups)) {
      barcode_groups <- apply(barcode_groups, 1, function(row) {
        barcode_group(
          name = row[["name"]],
          bc_list_name = row[["bc_list_name"]],
          max_edit_distance = ifelse("max_edit_distance" %in% names(row), as.numeric(row[["max_edit_distance"]]), 2)
        )
      })
    }
  }

  # Resolve bc_list_name keys/indices -> file paths (or fall back to barcodes_file)
  segments <- .resolve_bc_lists(
    segments,
    barcodes_files = barcodes_files
  )

  barcode_groups <- .resolve_bc_lists(
    barcode_groups,
    barcodes_files = barcodes_files
  )

  # Single sample
  stopifnot(length(fastq) == 1)
  if (file_test("-d", fastq)) {
    reads_in <- list.files(fastq, "\\.(fastq)|(fq)|(fasta)|(fa)$", full.names = TRUE)
  } else {
    reads_in <- fastq
  }
  read_counts <- flexiplex(
    r_segments = segments, r_barcode_groups = barcode_groups,
    max_flank_editdistance = max_flank_editdistance,
    reads_in = reads_in, 
    reads_out = reads_out, stats_out = stats_out, bc_out = tempfile(),
    reverseCompliment = reverseCompliment,
    n_threads = threads
  )

  # Cutadapt TSO trimming
  report <- list(read_counts = read_counts, stats_out = stats_out)
  if (is.character(TSO_seq) && nchar(TSO_seq) > 0 && TSO_prime %in% c(3, 5)) {
    # move the original reads to untrimmed_reads
    # run cutadapt and output to reads_out
    # move the untrimmed reads to noTSO_reads if full_length_only

    untrimmed_reads <- file.path(
      dirname(reads_out),
      paste0("untrimmed_", basename(reads_out))
    )
    noTSO_reads <- file.path(
      dirname(reads_out),
      paste0("noTSO_", basename(reads_out))
    )
    stopifnot(file.rename(reads_out, untrimmed_reads))
    tmp_json_file <- tempfile(fileext = ".json", tmpdir = dirname(reads_out))
    # cutadapt -a 'TSO' -o reads_out --untrimmed-output noTSO_out in_fq(untrimmed.fq)
    cutadapt(
      c(
        ifelse(TSO_prime == 3, "-a", "-g"),
        TSO_seq,
        "-o",
        reads_out,
        untrimmed_reads,
        "--json", tmp_json_file,
        if (full_length_only) {
          c("--untrimmed-output", noTSO_reads)
        },
        if (cutadapt_minimum_length > 0) {
          c("--minimum-length", cutadapt_minimum_length)
        }
      )
    )
    # parse cutadapt json report
    cutadapt_stats <- jsonlite::fromJSON(readLines(tmp_json_file, warn = FALSE))
    # cutadapt_stats <- cutadapt_stats[c('read1_with_adapter', 'input', 'output')]
    # cutadapt_stats <- unlist(cutadapt_stats)
    report$cutadapt <- cutadapt_stats
    unlink(tmp_json_file)
  } else {
    cat("Skipping TSO trimming...\n")
  }
  return(report)
}

#' @importFrom utils read.delim
convert_cellranger_bc <- function(bc_allow, bc_from, bc_to) {
  from <- read.delim(bc_from, header = FALSE)$V1
  to <- read.delim(bc_to, header = FALSE)$V1
  allowed <- read.delim(bc_allow, header = FALSE)$V1
  stopifnot(length(from) == length(to))
  stopifnot(all(allowed %in% from))
  return(to[from %in% allowed])
}

#' Plot Cell Barcode demultiplex statistics
#'
#' @description produce a barplot of cell barcode demultiplex statistics
#'
#' @importFrom dplyr bind_rows group_by summarise n_distinct ungroup arrange desc mutate row_number n
#' @importFrom ggplot2 ggplot aes geom_smooth geom_point scale_x_log10 scale_y_log10 theme_minimal ylab xlab ggtitle theme 
#' @param find_barcode_result output from \code{\link{find_barcode}}
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' fastq_dir <- tempfile()
#' dir.create(fastq_dir)
#' file.copy(system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
#'   file.path(fastq_dir, "musc_rps24.fastq.gz"))
#' sampled_lines <- readLines(file.path(fastq_dir, "musc_rps24.fastq.gz"), n = 400)
#' writeLines(sampled_lines, file.path(fastq_dir, "copy.fastq"))
#' bc_allow <- file.path(outdir, "bc_allow.tsv")
#' R.utils::gunzip(
#'   filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
#'   destname = bc_allow, remove = FALSE
#' )
#' find_barcode(
#'   fastq = fastq_dir,
#'   stats_out = file.path(outdir, "bc_stat.tsv.gz"),
#'   reads_out = file.path(outdir, "demultiplexed.fq"),
#'   barcodes_file = bc_allow, TSO_seq = "CCCATGTACTCTGCGTTGATACCACTGCTT"
#' ) |>
#'   plot_demultiplex_raw()
#' @return a list of ggplot objects:
#' \itemize{
#' \item reads_count_plot: stacked barplot of: demultiplexed reads
#' \item knee_plot: knee plot of UMI counts before TSO trimming
#' \item flank_editdistance_plot: flanking sequence (adaptor) edit-distance plot
#' \item barcode_editdistance_plot: barcode edit-distance plot
#' \item cutadapt_plot: if TSO trimming is performed, number of reads kept by cutadapt
#' }
#' @md
#' @export
#' @keywords internal
plot_demultiplex_raw <- function(find_barcode_result) {

  if (!all(sapply(find_barcode_result, is.list))) {
    find_barcode_result <- list("Sample" = find_barcode_result)
  }

  # stats file columns
  # readr's guessing will fail if the file is empty (no reads)
  stats_col_types <- list(
    "Read" = readr::col_character(),
    "CB" = readr::col_character(),
    "FlankEditDist" = readr::col_integer(),
    # "BarcodeEditDist" = readr::col_integer(),
    "UB" = readr::col_character()
    # "TooShort" = readr::col_logical()
  )

  find_barcode_result <- sapply(find_barcode_result, function(x) {
    if (is.null(x$reads_tb)) {
      x$reads_tb <- readr::read_delim(
        x$stats_out, delim = "\t", col_types = stats_col_types
      )
    }
    x
  }, simplify = FALSE)


  knee_tb <- sapply(find_barcode_result, "[[", "reads_tb", simplify = FALSE) |>
    dplyr::bind_rows(.id = "Sample") |>
    dplyr::mutate(Sample = factor(Sample)) |>
    dplyr::group_by(CB, Sample) |>
    dplyr::summarise(
      UMI_count = dplyr::n_distinct(UB)
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(Sample) |>
    dplyr::arrange(dplyr::desc(UMI_count)) |>
    dplyr::mutate(barcode_rank = dplyr::row_number()) |>
    dplyr::ungroup()


  # remove color if only one sample
  if (length(find_barcode_result) == 1) {
    knee_plot <- knee_tb |>
      ggplot2::ggplot(
        ggplot2::aes(x = barcode_rank, y = UMI_count)
      )
  } else {
    knee_plot <- knee_tb |>
      ggplot2::ggplot(
        ggplot2::aes(x = barcode_rank, y = UMI_count, col = Sample)
      )
  }

  knee_plot <- knee_plot +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_log10() +
    ggplot2::theme_minimal() +
    ggplot2::xlab("Barcode rank") +
    ggplot2::ylab("UMI counts (before TSO trimming)")
  if (nrow(knee_tb) > 10000) {
    if (requireNamespace("ggrastr", quietly = TRUE)) {
      knee_plot <- knee_plot +
        ggrastr::geom_point_rast(size = 0.5)
    } else {
      warning(
        "ggrastr package is not installed,",
        " using geom_point instead,",
        " this may be slow for large datasets",
        call. = FALSE
      )
      knee_plot <- knee_plot +
        ggplot2::geom_point(size = 0.5)
    }
  } else {
    knee_plot <- knee_plot +
      ggplot2::geom_smooth(
        span = 0.1, se = FALSE, method = "loess", alpha = 0.5, linewidth = 0.5
      ) +
      ggplot2::geom_point(size = 0.5)
  }

  flank_editdistance_plot <- sapply(
    find_barcode_result, "[[", "reads_tb", simplify = FALSE
  ) |>
    dplyr::bind_rows(.id = "Sample") |>
    dplyr::mutate(Sample = factor(Sample)) |>
    dplyr::group_by(FlankEditDist, Sample) |>
    dplyr::summarise(n_reads = dplyr::n()) |>
    dplyr::ungroup() |>
    ggplot2::ggplot(ggplot2::aes(x = FlankEditDist, y = n_reads, col = Sample)) +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() +
    ggplot2::ylab("number of reads") +
    ggplot2::xlab("Flanking sequence (adaptor) editdistance")
  if (length(find_barcode_result) == 1) {
    flank_editdistance_plot <- flank_editdistance_plot +
      ggplot2::theme(legend.position = "none")
  }

  # grouped (by sample) barplot
  # TODO: add back after implemented in flexiplex
  # barcode_editdistance_plot <- sapply(
  #   find_barcode_result, "[[", "reads_tb", simplify = FALSE
  # ) |>
  #   dplyr::bind_rows(.id = "Sample") |>
  #   dplyr::mutate(Sample = factor(Sample)) |>
  #   dplyr::group_by(BarcodeEditDist, Sample) |>
  #   dplyr::ungroup() |>
  #   dplyr::mutate(BarcodeEditDist = factor(BarcodeEditDist)) |>
  #   ggplot2::ggplot(ggplot2::aes(x = Sample, fill = BarcodeEditDist)) +
  #   ggplot2::geom_bar(stat = "count", position = "dodge") +
  #   ggplot2::theme_minimal() +
  #   ggplot2::ylab("number of reads")
  # if (length(find_barcode_result) == 1) {
  #   barcode_editdistance_plot <- barcode_editdistance_plot +
  #     ggplot2::theme(legend.position = "none")
  # }

  # read counts
  read_counts_plot <-
    sapply(find_barcode_result, function(x) x$read_counts, simplify = FALSE) |>
      dplyr::bind_rows(.id = "Sample") |>
      dplyr::mutate(`undemultiplexted reads` = `total reads` - `demultiplexed reads`,
        `multi-matching reads` = `demultiplexed reads` - `single match reads`,
      ) |>
      dplyr::select(Sample, `undemultiplexted reads`, `multi-matching reads`, `single match reads`) |>
      tidyr::pivot_longer(
        !Sample,
        names_to = "Type", values_to = "Reads"
      ) |> # make stacked barplot, each sample is a bar
      dplyr::mutate(Type = factor(Type, levels = c(
        "undemultiplexted reads",
        "multi-matching reads",
        "single match reads"
      ))) |>
      ggplot2::ggplot(ggplot2::aes(x = Sample, y = Reads, fill = Type)) +
      ggplot2::geom_col(position = "stack") +
      ggplot2::scale_fill_manual(values = c(
        "undemultiplexted reads" = "#B3B3B3",
        "single match reads" = "#8DA0CB",
        "multi-matching reads" = "#FC8D62")) +
      ggplot2::theme_minimal() +
      ggplot2::ylab("number of reads") +
      ggplot2::xlab("Sample") +
      ggplot2::ggtitle("Read counts")

  # Cutadapt report
  if (all(sapply(find_barcode_result, function(x) "cutadapt" %in% names(x)))) {
    cutadapt_plot <- 
      sapply(find_barcode_result, function(x) x$cutadapt$read_counts, simplify = FALSE) |>
      dplyr::bind_rows(.id = "Sample") |>
      tidyr::pivot_longer(
        cols = c(input, output, read1_with_adapter),
        names_to = "Type", values_to = "Count"
      ) |>
      ggplot2::ggplot(ggplot2::aes(x = Sample, y = Count, fill = Type)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::theme_minimal() +
      ggplot2::ylab("number of reads") +
      ggplot2::xlab("Sample") +
      ggplot2::ggtitle("Cutadapt results")
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  return(list(
    reads_count_plot = read_counts_plot,
    knee_plot = knee_plot,
    flank_editdistance_plot = flank_editdistance_plot,
    # barcode_editdistance_plot = barcode_editdistance_plot,
    cutadapt_plot = switch(exists("cutadapt_plot"), cutadapt_plot, NULL)
  ))
}

default_names <- function(x, name) {
  if (is.null(names(x))) {
    names(x) <- name
  }
  return(x)
}
default_values <- function(x, value) {
  if (length(x) != length(value)) {
    x <- value
  }
  return(x)
}
