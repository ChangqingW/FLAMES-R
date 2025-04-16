#' Create Configuration File From Arguments
#'
#' @details Create a list object containing the arguments supplied in a format usable for the FLAMES pipeline.
#' Also writes the object to a JSON file, which is located with the prefix 'config_' in the supplied \code{outdir}.
#' Default values from \code{extdata/config_sclr_nanopore_3end.json} will be used for unprovided parameters.
#'
#' @param outdir the destination directory for the configuratio nfile
#' @param type use an example config, available values:
#' \describe{
#'  \item{"sc_3end"}{ - config for 10x 3' end ONT reads}
#'  \item{"SIRV"}{ - config for the SIRV example reads}
#' }
#' @param ... Configuration parameters.
#' \describe{
#'  \item{seed}{ - Integer. Seed for minimap2.}
#'  \item{threads}{ - Number of threads to use.}
#'  \item{do_barcode_demultiplex}{ - Boolean. Specifies whether to run the barcode demultiplexing step.}
#'  \item{do_genome_alignment}{ - Boolean. Specifies whether to run the genome alignment step. \code{TRUE} is recommended}
#'  \item{do_gene_quantification}{ - Boolean. Specifies whether to run gene quantification using the genome alignment results. \code{TRUE} is recommended}
#'  \item{do_isoform_identification}{ - Boolean. Specifies whether to run the isoform identification step. \code{TRUE} is recommended}
#'  \item{bambu_isoform_identification}{ - Boolean. Whether to use Bambu for isoform identification.}
#'  \item{multithread_isoform_identification}{ - Boolean. Whether to use FLAMES' new multithreaded Cpp implementation for isoform identification.}
#'  \item{do_read_realignment}{ - Boolean. Specifies whether to run the read realignment step. \code{TRUE} is recommended}
#'  \item{do_transcript_quantification}{ - Boolean. Specifies whether to run the transcript quantification step. \code{TRUE} is recommended}
#'  \item{barcode_parameters}{ - List. Parameters for barcode demultiplexing passed to \code{find_barcode} (except \code{fastq}, \code{barcodes_file}, \code{stats_out}, \code{reads_out}) and \code{threads}, which are set by the pipeline, see \code{?find_barcode} for more details.}
#'  \item{generate_raw_isoform}{ - Boolean. Whether to generate all isoforms for debugging purpose.}
#'  \item{max_dist}{ - Maximum distance allowed when merging splicing sites in isoform consensus clustering.}
#'  \item{max_ts_dist}{ - Maximum distance allowed when merging transcript start/end position in isoform consensus clustering.}
#'  \item{max_splice_match_dist}{ - Maximum distance allowed when merging splice site called from the data and the reference annotation.}
#'  \item{min_fl_exon_len}{ - Minimum length for the first exon outside the gene body in reference annotation. This is to correct the alignment artifact}
#'  \item{max_site_per_splice}{ - Maximum transcript start/end site combinations allowed per splice chain}
#'  \item{min_sup_cnt}{ - Minimum number of read support an isoform decrease this number will significantly increase the number of isoform detected.}
#'  \item{min_cnt_pct}{ - Minimum percentage of count for an isoform relative to total count for the same gene.}
#'  \item{min_sup_pct}{ - Minimum percentage of count for an splice chain that support a given transcript start/end site combination.}
#'  \item{strand_specific}{ - 0, 1 or -1. 1 indicates if reads are in the same strand as mRNA, -1 indicates reads are reverse complemented, 0 indicates reads are not strand specific.}
#'  \item{remove_incomp_reads}{ - The strenge of truncated isoform filtering. larger number means more stringent filtering.}
#'  \item{use_junctions}{ - whether to use known splice junctions to help correct the alignment results}
#'  \item{no_flank}{ - Boolean. for synthetic spike-in data. refer to Minimap2 document for detail}
#'  \item{use_annotation}{ - Boolean. whether to use reference to help annotate known isoforms}
#'  \item{min_tr_coverage}{ - Minimum percentage of isoform coverage for a read to be aligned to that isoform}
#'  \item{min_read_coverage}{ - Minimum percentage of read coverage for a read to be uniquely aligned to that isoform}
#' }
#'
#' @return file path to the config file created
#' @examples
#' # create the default configuration file
#' outdir <- tempdir()
#' config <- create_config(outdir)
#'
#' @importFrom jsonlite toJSON fromJSON
#' @importFrom utils modifyList
#' @importFrom stats setNames
#' @export
# Centralized default configuration
get_default_config <- function() {
  jsonlite::fromJSON(
    system.file("extdata", "config_sclr_nanopore_3end.json", package = "FLAMES")
  )
}

# Utility function to recursively merge nested lists
merge_nested_lists <- function(default, updates) {
  for (name in names(updates)) {
    if (is.list(updates[[name]]) && is.list(default[[name]])) {
      default[[name]] <- merge_nested_lists(default[[name]], updates[[name]])
    } else {
      default[[name]] <- updates[[name]]
    }
  }
  return(default)
}

# Updated create_config function
create_config <- function(outdir, type = "sc_3end", ...) {
  # Load default configuration
  config <- get_default_config()

  # Apply type-specific modifications
  if (type == "SIRV") {
    config$alignment_parameters$no_flank <- TRUE
  } else if (type != "sc_3end") {
    stop("Unrecognised config type ", type)
  }

  # Merge user-specified updates
  updates <- list(...)
  if (length(updates) > 0) {
    if (any(is.null(names(updates))) || "" %in% names(updates)) {
      stop("Parameters must be named")
    }
    config <- merge_nested_lists(config, updates)
  }

  # Write created config file
  config_file_path <- file.path(outdir, paste0("config_file_", Sys.getpid(), ".json"))
  cat("Writing configuration parameters to: ", config_file_path, "\n")
  write(jsonlite::toJSON(config, pretty = TRUE), config_file_path)

  return(config_file_path)
}

# Define the base S4 class for FLAMES pipelines
setClass(
  "FLAMES.Pipeline",
  slots = list(
    config = "list",             # Configuration parameters
    steps = "character",        # Steps to perform
    completed_steps = "character", # Completed steps
    metadata = "list",          # Metadata for the pipeline run
    outdir = "character"        # Output directory
  )
)

# Define a method to initialize the pipeline
setGeneric("initializePipeline", function(object, config, steps, outdir) {
  standardGeneric("initializePipeline")
})

setMethod(
  "initializePipeline",
  "FLAMES.Pipeline",
  function(object, config, steps, outdir) {
    object@config <- config
    object@steps <- steps
    object@completed_steps <- character(0)
    object@metadata <- list()
    object@outdir <- outdir
    return(object)
  }
)

#' @importFrom Matrix tail
#' @importFrom stringr str_split
#' @importFrom jsonlite fromJSON
# Updated check_arguments function to ensure backward compatibility
check_arguments <- function(
    annotation, fastq, genome_bam,
    outdir, genome_fa, config_file) {
  if (!dir.exists(outdir)) {
    cat("Output directory does not exists: one is being created\n")
    dir.create(outdir)
    print(outdir)
  }

  if (is.null(config_file)) {
    cat("No config file provided, creating a default config in", outdir, "\n")
    config_file <- create_config(outdir)
  }

  # Load configuration and ensure backward compatibility
  config <- jsonlite::fromJSON(config_file)
  default_config <- get_default_config()
  config <- merge_nested_lists(default_config, config)

  # Validate configuration
  if (config$isoform_parameters$downsample_ratio > 1 || config$isoform_parameters$downsample_ratio <= 0) {
    stop("downsample_ratio should be between 0 and 1")
  }
  if (!is.null(fastq) && any(!file.exists(fastq))) {
    stop(paste0("Make sure ", fastq, " exists."))
  }
  if (!file.exists(annotation)) {
    stop(paste0("Make sure ", annotation, " exists."))
  }
  if (!file.exists(genome_fa)) {
    stop(paste0("Make sure ", genome_fa, " exists."))
  }
  if (!is.null(genome_bam) && any(!file.exists(genome_bam))) {
    stop("Make sure genome_bam exists")
  }

  if (config$pipeline_parameters$bambu_isoform_identification) {
    if (!(stringr::str_ends(annotation, ".gtf") | stringr::str_ends(annotation, ".gtf.gz"))) {
      stop("Bambu requires GTF format for annotation file.\n")
    }
  }

  if (config$pipeline_parameters$do_transcript_quantification &&
    config$pipeline_parameters$oarfish_quantification &&
    !config$pipeline_parameters$do_gene_quantification) {
    warning("You have set to use oarfish quantification without gene quantification. Oarfish currently does not collapse UMIs, and gene quantification performs UMI collapsing. You may want to set do_gene_quantification to TRUE for more accurate results.")
  }

  return(list(config = config))
}
