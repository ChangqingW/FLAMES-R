#' Convert GFF/GTF to BED file
#'
#' @param gff Path to the GFF/GTF file
#' @param bed Path to the output BED file to be written
#'
#' @return invisible, the BED file is written to the specified path
#' @importFrom rtracklayer import export.bed
#' @keywords internal
gff2bed <- function(gff, bed) {
  gff <- rtracklayer::import(gff)
  gff$score <- as.numeric(rep(0, length(gff)))
  return(invisible(rtracklayer::export.bed(gff, bed)))
}

#' Minimap2 Align to Genome
#'
#' @description
#' Uses minimap2 to align sequences agains a reference databse.
#' Uses options '-ax splice -t 12 -k14 --secondary=no \code{fa_file} \code{fq_in}'
#'
#' @param fq_in File path to the fastq file used as a query sequence file
#' @param fa_file Path to the fasta file used as a reference database for alignment
#' @param config Parsed list of FLAMES config file
#' @param outfile Path to the output file
#' @param minimap2_args Arguments to pass to minimap2, see minimap2 documentation for details.
#' @param sort_by Column to sort the bam file by, see \code{samtools sort} for details
#' @param minimap2 Path to minimap2 binary
#' @param samtools path to the samtools binary.
#' @param threads Integer, threads for minimap2 to use, see minimap2 documentation for details,
#' @param tmpdir Temporary directory to use for intermediate files.
#' FLAMES will try to detect cores if this parameter is not provided.
#'
#' @return a \code{data.frame} summarising the reads aligned
#'
#' @importFrom Rsamtools sortBam indexBam asBam
#' @keywords internal
minimap2_align <- function(fq_in, fa_file, config, outfile, minimap2_args, sort_by, minimap2, samtools, threads = 1, tmpdir) {
  has_tags <- grepl("\t", readLines(fq_in, n = 1))
  if (has_tags && !any(grepl("-y", minimap2_args))) {
    message("Your fastq file appears to have tags, but you did not provide the -y option to minimap2 to include the tags in the output.")
  }

  tmp_bam <- tempfile(fileext = ".bam", tmpdir = tmpdir)
  if (is.character(samtools) && !is.na(samtools) && file.exists(samtools)) {
    have_samtools <- TRUE
    minimap2_status <- base::system2(
      command = minimap2,
      args = base::append(
        minimap2_args,
        c(fa_file, fq_in, "|", samtools, "view -b -o", tmp_bam, "-")
      )
    )
    if (!is.null(base::attr(minimap2_status, "status")) &&
      base::attr(minimap2_status, "status") != 0) {
      stop(paste0("error running minimap2:\n", minimap2_status))
    }
  } else {
    warning("samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.")
    have_samtools <- FALSE
    tmp_sam <- tempfile(fileext = ".sam", tmpdir = tmpdir)
    minimap2_status <- base::system2(
      command = minimap2,
      args = base::append(
        minimap2_args,
        c(fa_file, fq_in, ">", tmp_sam)
      )
    )
    if (!is.null(base::attr(minimap2_status, "status")) &&
      base::attr(minimap2_status, "status") != 0) {
      stop(paste0("error running minimap2:\n", minimap2_status))
    }
    Rsamtools::asBam(
      tmp_sam,
      # great, destination is not the destination in Rsamtools
      destination = sub("\\.bam", "", tmp_bam),
      overwrite = TRUE,
      indexDestination = FALSE
    )
    unlink(tmp_sam)
  }

  tmp_bam_sorted <- tempfile(fileext = ".sort", tmpdir = tmpdir)
  sort_args <- c(
    "sort",
    "-@", threads,
    "-T", tmp_bam_sorted,
    tmp_bam, "-o", outfile
  )

  if (sort_by != "none") {
    if (sort_by == "coordinates") {
      message(sprintf("Sorting BAM files by genome coordinates with %s threads...\n", threads))
    } else {
      sort_args <- c(sort_args, "-t", sort_by)
      message(sprintf("Sorting BAM files by %s with %s threads...\n", threads, sort_by))
    }

    if (have_samtools) {
      sort_status <- base::system2(
        command = samtools,
        args = sort_args
      )
      if (!is.null(base::attr(sort_status, "status")) &&
        base::attr(sort_status, "status") != 0) {
        stop(paste0("error running samtools sort:\n", sort_status))
      }
    } else {
      Rsamtools::sortBam(
        file = tmp_bam,
        destination = sub("\\.bam", "", outfile),
        byTag = if (sort_by %in% c("coordinates", "none")) NULL else sort_by,
        byQname = sort_by == "none",
        nThreads = threads
      )
    }
  } else {
    message("Skipped sorting BAM files.\n")
    # most platforms will not rename files
    # from one file system to another
    # and file.rename will fail
    # file.xxx functions fail witout throwing errors
    # and return FALSE
    if (!file.rename(tmp_bam, outfile)) {
      file.copy(tmp_bam, outfile)
    }
  }

  if (sort_by == "coordinates") {
    cat("Indexing bam files\n")
    if (have_samtools) {
      index_status <- base::system2(
        command = samtools,
        args = c("index", outfile)
      )
      if (!is.null(base::attr(index_status, "status")) &&
        base::attr(index_status, "status") != 0) {
        stop(paste0("error running samtools index:\n", index_status))
      }
    } else {
      Rsamtools::indexBam(outfile)
    }
  }

  unlink(c(tmp_bam, tmp_bam_sorted))
  return(get_flagstat(outfile, samtools))
}

#' @importFrom Rsamtools quickBamFlagSummary
get_flagstat <- function(bam, samtools_path) {
  stats_df <- data.frame(total = 0, mapped = 0, primary = 0, secondary = 0)
  rownames(stats_df) <- bam
  if (!missing("samtools_path") && !is.na(samtools_path) && file.exists(samtools_path)) {
    output <- base::system2(samtools_path, c("flagstat", bam), stdout = TRUE)
    stats_df["total"] <- as.numeric(regmatches(
      output[grepl("in total", output)],
      regexpr("\\d+", output[grepl("in total", output)])
    )[1])
    stats_df["mapped"] <- as.numeric(regmatches(
      output[grepl("mapped", output)],
      regexpr("(?<!S)\\d+", output[grepl("mapped", output)], perl = TRUE)
    )[1])
    stats_df["primary"] <- as.numeric(regmatches(
      output[grepl("primary$", output)],
      regexpr("(?<!S)\\d+", output[grepl("primary$", output)], perl = TRUE)
    )[1])
    stats_df["secondary"] <- as.numeric(regmatches(output[grepl(
      "secondary$",
      output
    )], regexpr("(?<!S)\\d+", output[grepl("secondary$", output)],
      perl = TRUE
    ))[1])
  } else {
    output <- utils::capture.output(Rsamtools::quickBamFlagSummary(bam))
    stats_df["total"] <- as.numeric(regmatches(
      output[grepl("^All records", output)],
      regexpr("\\d+", output[grepl("^All records", output)])
    )[1])
    stats_df["mapped"] <- as.numeric(regmatches(output[grepl(
      "record is mapped",
      output
    )], regexpr("(?<!S)\\d+", output[grepl("record is mapped", output)],
      perl = TRUE
    ))[1])
    stats_df["primary"] <- as.numeric(regmatches(output[grepl(
      "primary alignment",
      output
    )], regexpr("(?<!S)\\d+", output[grepl("primary alignment", output)],
      perl = TRUE
    ))[1])
    stats_df["secondary"] <- as.numeric(regmatches(output[grepl(
      "secondary alignment",
      output
    )], regexpr("(?<!S)\\d+", output[grepl("primary alignment", output)],
      perl = TRUE
    ))[1])
  }
  stats_df
}

# generic pie chart
#' @importFrom ggplot2 ggplot aes geom_bar ggtitle coord_polar element_blank theme position_stack theme_bw geom_histogram ggtitle ylab xlab geom_text
plot_flagstat <- function(flagstat) {
  flagstat[["unmapped"]] <- flagstat[["total"]] - flagstat[["mapped"]]
  tidyr::pivot_longer(as.data.frame(flagstat), everything()) |>
    dplyr::filter(!name %in% c("total", "mapped")) |>
    ggplot(aes(x = "", y = value, label = value, fill = name)) +
    geom_bar(stat = "identity") +
    coord_polar("y") +
    ggtitle("Alignment summary") +
    geom_text(position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}
