#' Minimap2 Align to Genome
#'
#' @description
#' Uses minimap2 to align sequences agains a reference databse.
#' Uses options '-ax splice -t 12 -k14 --secondary=no \code{fa_file} \code{fq_in}'
#'
#' @param config Parsed list of FLAMES config file
#' @param fa_file Path to the fasta file used as a reference database for alignment
#' @param fq_in File path to the fastq file used as a query sequence file
#' @param outdir Output folder
#' @param minimap2 Path to minimap2 binary
#' @param k8 Path to the k8 Javascript shell binary
#' @param prefix String, the prefix (e.g. sample name) for the outputted BAM file
#' @param threads Integer, threads for minimap2 to use, see minimap2 documentation for details,
#' FLAMES will try to detect cores if this parameter is not provided.
#' @param samtools path to the samtools binary, required for large datasets since \code{Rsamtools} does not support \code{CSI} indexing
#'
#' @return a \code{data.frame} summarising the reads aligned
#' @seealso [minimap2_realign()]
#'
#' @importFrom Rsamtools sortBam indexBam asBam
#' @export
#' @examples
#' temp_path <- tempfile()
#' bfc <- BiocFileCache::BiocFileCache(temp_path, ask = FALSE)
#' file_url <- 'https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data'
#' fastq1 <- bfc[[names(BiocFileCache::bfcadd(bfc, 'Fastq1', paste(file_url, 'fastq/sample1.fastq.gz', sep = '/')))]]
#' genome_fa <- bfc[[names(BiocFileCache::bfcadd(bfc, 'genome.fa', paste(file_url, 'SIRV_isoforms_multi-fasta_170612a.fasta', sep = '/')))]]
#' annotation <- bfc[[names(BiocFileCache::bfcadd(bfc, 'annot.gtf', paste(file_url, 'SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf', sep = '/')))]]
#' outdir <- tempfile()
#' dir.create(outdir)
#' minimap2_align(
#'   config = jsonlite::fromJSON(
#'     system.file("extdata", "config_sclr_nanopore_3end.json", package = 'FLAMES')
#'   ),
#'   fa_file = genome_fa,
#'   fq_in = fastq1,
#'   outdir = outdir
#' )
minimap2_align <- function(fq_in, fa_file, config, outfile, minimap2_args, sort_by, minimap2, samtools, threads = 1) {

  has_tags <- grepl("\t", readLines(fq_in, n = 1))
  if (has_tags && !any(grepl("-y", minimap2_args))) {
    message("Your fastq file appears to have tags, but you did not provide the -y option to minimap2 to include the tags in the output.")
  }

  # TODO: use Rsamtools if samtools is not available
  stopifnot("Samtools not found" = !is.na(samtools))
  tmp_bam <- tempfile(fileext = ".bam")
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

  sort_args <- c(
    "sort",
    "-@", threads,
    "-T", tempfile(fileext = ".sort"),
    tmp_bam, "-o", outfile
  )

  if (sort_by != "none") {
    if (sort_by == "coordinates") {
      message(sprintf("Sorting BAM files by genome coordinates with %s threads...\n", threads))
    } else {
      sort_args <- c(sort_args, "-t", sort_by)
      message(sprintf("Sorting BAM files by %s with %s threads...\n", threads, sort_by))
    }

    sort_status <- base::system2(
      command = samtools,
      args = sort_args
    )
    if (!is.null(base::attr(sort_status, "status")) &&
          base::attr(sort_status, "status") != 0) {
      stop(paste0("error running samtools sort:\n", sort_status))
    }

  } else {
    message("Skipped sorting BAM files.\n")
    file.rename(tmp_bam, outfile)
  }

  if (sort_by == "coordinates") {
    cat("Indexing bam files\n")
    index_status <- base::system2(
      command = samtools,
      args = c("index", outfile)
    )
    if (!is.null(base::attr(index_status, "status")) &&
          base::attr(index_status, "status") != 0) {
      stop(paste0("error running samtools index:\n", index_status))
    }
  }

  if (!is.na(samtools)) {
    return(get_flagstat(outfile, samtools))
  }
  # No equivalent to samtools flagstat in Rsamtools
  # Rsamtools::quickBamFlagSummary does not return anything
}


#' Find path to a binary
#' Wrapper for Sys.which to find path to a binary
#' @importFrom withr with_path
#' @importFrom basilisk obtainEnvironmentPath
#' @description
#' This function is a wrapper for \code{base::Sys.which} to find the path
#' to a command. It also searches within the \code{FLAMES} basilisk conda
#' environment. This function also replaces "" with \code{NA} in the
#' output of \code{base::Sys.which} to make it easier to check if the
#' binary is found.
#' @param command character, the command to search for
#' @return character, the path to the command or \code{NA}
#' @examples
#' find_bin("minimap2")
#' @export
find_bin <- function(command) {
  conda_bins <- file.path(basilisk::obtainEnvironmentPath(bins_env), 'bin')
  which_command <- withr::with_path(
    new = conda_bins,
    action = "suffix",
    code = Sys.which(command)
  )
  # replace "" with NA
  which_command[which_command == ""] <- NA
  return(which_command)
}

# total mapped primary secondary
get_flagstat <- function(bam, samtools_path) {
  stats_df <- data.frame(total = 0, mapped = 0, primary = 0, secondary = 0)
  rownames(stats_df) <- bam
  if (!missing("samtools_path") && is.character(samtools_path)) {
    output <- base::system2(samtools_path, c("flagstat", bam), stdout = TRUE)
    stats_df["total"] <- as.numeric(regmatches(output[grepl("in total", output)],
      regexpr("\\d+", output[grepl("in total", output)]))[1])
    stats_df["mapped"] <- as.numeric(regmatches(output[grepl("mapped", output)],
      regexpr("(?<!S)\\d+", output[grepl("mapped", output)], perl = TRUE))[1])
    stats_df["primary"] <- as.numeric(regmatches(output[grepl("primary$", output)],
      regexpr("(?<!S)\\d+", output[grepl("primary$", output)], perl = TRUE))[1])
    stats_df["secondary"] <- as.numeric(regmatches(output[grepl("secondary$",
      output)], regexpr("(?<!S)\\d+", output[grepl("secondary$", output)],
      perl = TRUE))[1])
  } else {
    output <- utils::capture.output(Rsamtools::quickBamFlagSummary(bam))
    stats_df["total"] <- as.numeric(regmatches(output[grepl("^All records", output)],
      regexpr("\\d+", output[grepl("^All records", output)]))[1])
    stats_df["mapped"] <- as.numeric(regmatches(output[grepl("record is mapped",
      output)], regexpr("(?<!S)\\d+", output[grepl("record is mapped", output)],
      perl = TRUE))[1])
    stats_df["primary"] <- as.numeric(regmatches(output[grepl("primary alignment",
      output)], regexpr("(?<!S)\\d+", output[grepl("primary alignment", output)],
      perl = TRUE))[1])
    stats_df["secondary"] <- as.numeric(regmatches(output[grepl("secondary alignment",
      output)], regexpr("(?<!S)\\d+", output[grepl("primary alignment", output)],
      perl = TRUE))[1])
  }
  stats_df
}

# generic pie chart
#' @importFrom ggplot2 ggplot aes geom_bar ggtitle coord_polar element_blank theme position_stack theme_bw geom_histogram ggtitle ylab xlab geom_text
plot_flagstat <- function(flagstat) {
  flagstat[["unmapped"]] <- flagstat[["total"]] - flagstat[["mapped"]]
  tidyr::pivot_longer(as.data.frame(flagstat), everything()) %>%
    dplyr::filter(!name %in% c("total", "mapped")) %>%
    ggplot(aes(x = "", y = value, label = value, fill = name)) + geom_bar(stat = "identity") +
    coord_polar("y") + ggtitle("Alignment summary") + geom_text(position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL) + theme_bw() + theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(),
      axis.ticks.y = element_blank())
}
