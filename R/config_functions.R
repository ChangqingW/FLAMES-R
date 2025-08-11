#' Load Configurations
#'
#' Loads a configuration file and fills in missing values with defaults
#' from the package's default configuration.
#'
#' @param config_file Path to the configuration JSON file
#' @param type Config type to use for defaults ("sc_3end" or "SIRV")
#'
#' @return A complete configuration list with all parameters filled
#' @importFrom jsonlite fromJSON
#' @importFrom utils modifyList
#' @export
load_config <- function(config_file, type = "sc_3end") {
  # Load default configuration
  if (type == "sc_3end") {
    default_config <- jsonlite::fromJSON(
      system.file(
        "extdata", "config_sclr_nanopore_3end.json", package = "FLAMES"
      )
    )
  } else if (type == "SIRV") {
    default_config <- jsonlite::fromJSON(
      system.file(
        "extdata", "config_sclr_nanopore_3end.json", package = "FLAMES"
      )
    )
    default_config$alignment_parameters$no_flank <- TRUE
  } else {
    stop("Unrecognised config type ", type)
  }

  # Load user configuration
  user_config <- jsonlite::fromJSON(config_file)

  # Recursively merge configurations, with user config taking precedence
  # but missing values filled from defaults
  merged_config <- merge_configs_recursive(default_config, user_config)

  return(merged_config)
}

#' Recursively Merge Configuration Lists
#'
#' Internal function to recursively merge configuration lists,
#' filling missing values from defaults while preserving user values
#'
#' @param default_config Default configuration list
#' @param user_config User configuration list
#' @return Merged configuration list
#' @note Special case: when user_config contains barcode_parameters.pattern
#'   as a list, the entire pattern list is preserved as-is without merging
#'   with defaults to maintain user-specified order and structure.
#' @keywords internal
merge_configs_recursive <- function(default_config, user_config) {
  # Remove comment field from default if it exists
  if ("comment" %in% names(default_config)) {
    default_config$comment <- NULL
  }

  # If user_config is NULL or empty, return default
  if (is.null(user_config) || length(user_config) == 0) {
    return(default_config)
  }

  # For each element in default_config
  for (name in names(default_config)) {
    if (name %in% names(user_config)) {
      # Special case: if we're dealing with barcode_parameters and user has pattern as a list,
      # don't merge the pattern, use it as-is
      if (name == "barcode_parameters" &&
          is.list(default_config[[name]]) &&
          is.list(user_config[[name]]) &&
          "pattern" %in% names(user_config[[name]]) &&
          is.list(user_config[[name]]$pattern)) {

        # Merge barcode_parameters but keep user pattern as-is
        merged_barcode_params <- merge_configs_recursive(
          default_config[[name]],
          user_config[[name]]
        )
        # Override the pattern with user's exact pattern (no merging)
        merged_barcode_params$pattern <- user_config[[name]]$pattern
        default_config[[name]] <- merged_barcode_params

      } else if (is.list(default_config[[name]]) && is.list(user_config[[name]])) {
        # If both are lists, merge recursively (normal case)
        default_config[[name]] <- merge_configs_recursive(
          default_config[[name]],
          user_config[[name]]
        )
      } else {
        # Use user value (non-list or user overrides)
        default_config[[name]] <- user_config[[name]]
      }
    }
    # If name not in user_config, keep default value (already there)
  }

  # Add any additional parameters from user_config that aren't in default
  for (name in names(user_config)) {
    if (!name %in% names(default_config)) {
      default_config[[name]] <- user_config[[name]]
    }
  }

  default_config
}

#' Set Nested Configuration Parameter
#'
#' Helper function to set a nested parameter in a configuration list
#' using dot notation (e.g., "barcode_parameters.pattern.primer")
#'
#' @param config Configuration list
#' @param param_path Parameter path using dot notation
#' @param value Value to set
#' @return Modified configuration list
#' @keywords internal
set_nested_param <- function(config, param_path, value) {
  path_parts <- strsplit(param_path, "\\.")[[1]]

  # Recursive helper function
  set_nested_recursive <- function(obj, parts, val) {
    if (length(parts) == 1) {
      obj[[parts[1]]] <- val
      return(obj)
    }

    first_part <- parts[1]
    remaining_parts <- parts[-1]

    # Initialize as empty list if doesn't exist
    if (!first_part %in% names(obj)) {
      obj[[first_part]] <- list()
    }

    # Recurse
    obj[[first_part]] <- set_nested_recursive(
      obj[[first_part]], remaining_parts, val
    )
    return(obj)
  }

  return(set_nested_recursive(config, path_parts, value))
}

#' Create Configuration File From Arguments
#'
#' @details Create a list object containing the arguments supplied in a format
#' usable for the FLAMES pipeline, and writes the object to a JSON file, which
#' is located with the prefix 'config_' in the supplied \code{outdir}. Default
#' values from \code{extdata/config_sclr_nanopore_3end.json} will be used for
#' unprovided parameters.
#'
#' Parameters can be specified using dot to indicate nested sections, e.g.,
#' \code{barcode_parameters.max_bc_editdistance = 3} or
#' \code{barcode_parameters.pattern.primer = "ATCG"}. Alternatively, you can
#' open the created config file and edit it manually.
#'
#' @param outdir the destination directory for the configuration file
#' @param type use an example config, available values:
#' \describe{
#'  \item{"sc_3end"}{ - config for 10x 3' end ONT reads}
#'  \item{"SIRV"}{ - config for the SIRV example reads}
#' }
#' @param ... Configuration parameters (using dot for nested parameters)
#' \describe{
#'  \item{seed}{ - Integer. Seed for minimap2.}
#'  \item{threads}{ - Number of threads to use.}
#'  \item{do_barcode_demultiplex}{ - Boolean. Specifies whether to run the
#'  barcode demultiplexing step.}
#'  \item{do_genome_alignment}{ - Boolean. Specifies whether to run the genome
#'  alignment step. \code{TRUE} is recommended}
#'  \item{do_gene_quantification}{ - Boolean. Specifies whether to run gene
#'  quantification using the genome alignment results. \code{TRUE} is
#'  recommended}
#'  \item{do_isoform_identification}{ - Boolean. Specifies whether to run the
#'  isoform identification step. \code{TRUE} is recommended}
#'  \item{bambu_isoform_identification}{ - Boolean. Whether to use Bambu for
#'  isoform identification.}
#'  \item{multithread_isoform_identification}{ - Boolean. Whether to use FLAMES'
#'  new multithreaded Cpp implementation for isoform identification.}
#'  \item{do_read_realignment}{ - Boolean. Specifies whether to run the read
#'  realignment step. \code{TRUE} is recommended}
#'  \item{do_transcript_quantification}{ - Boolean. Specifies whether to run the
#'  transcript quantification step. \code{TRUE} is recommended}
#'  \item{barcode_parameters.max_bc_editdistance}{ - Maximum edit distance for
#'  barcode matching}
#'  \item{barcode_parameters.pattern.primer}{ - Primer sequence pattern}
#'  \item{isoform_parameters.max_dist}{ - Maximum distance allowed when merging
#'  splicing sites}
#'  \item{...}{ - Other nested parameters, using dot to indicate nested section}
#' }
#'
#' @return file path to the config file created
#' @examples
#' # create the default configuration file
#' outdir <- tempdir()
#' config <- create_config(outdir)
#'
#' # create config with custom parameters including nested ones
#' config <- create_config(outdir,
#'   threads = 16,
#'   barcode_parameters.max_bc_editdistance = 3,
#'   barcode_parameters.pattern.primer = "ATCGATCG",
#'   isoform_parameters.min_sup_cnt = 10
#' )
#'
#' @importFrom jsonlite toJSON fromJSON
#' @export
create_config <- function(outdir, type = "sc_3end", ...) {
  # Load default configuration
  if (type == "sc_3end") {
    config <- jsonlite::fromJSON(system.file(
      "extdata", "config_sclr_nanopore_3end.json", package = "FLAMES"
    ))
  } else if (type == "SIRV") {
    config <- jsonlite::fromJSON(system.file(
      "extdata", "config_sclr_nanopore_3end.json", package = "FLAMES"
    ))
    config$alignment_parameters$no_flank <- TRUE
  } else {
    stop("Unrecognised config type ", type)
  }

  # Get user parameters
  updates <- list(...)

  if (length(updates) > 0) {
    if (any(is.null(names(updates))) || "" %in% names(updates)) {
      stop("Parameters must be named")
    }

    # Remove comment field
    config <- within(config, rm(comment))

    # Apply each parameter update
    for (param_name in names(updates)) {
      param_value <- updates[[param_name]]

      # Check if parameter uses dot notation for nested access
      if (grepl("\\.", param_name)) {
        config <- set_nested_param(config, param_name, param_value)
      } else {
        # Handle top-level parameters (maintain backward compatibility)
        # Try to find which section this parameter belongs to
        found <- FALSE
        for (section_name in names(config)) {
          if (is.list(config[[section_name]]) &&
                param_name %in% names(config[[section_name]])) {
            config[[section_name]][[param_name]] <- param_value
            found <- TRUE
            break
          }
        }

        # If not found in any section, provide helpful error with example
        if (!found) {
          # Show example using actual nested parameters from the config
          stop(
            "Unknown parameter '", param_name,
            "'. Nested Parameters should be specified",
            " using dots to indicate sections. ",
            "Examples:\n",
            "  - For 'threads' use: pipeline_parameters.threads\n",
            "  - For 'max_bc_editdistance' use:",
            " barcode_parameters.max_bc_editdistance\n",
            "Alternatively, you can open the created config file ",
            "and edit it manually.\n",
          )
        }
      }
    }
  }

  # Write created config file
  config_file_path <- file.path(
    outdir, paste0("config_file_", Sys.getpid(), ".json")
  )
  cat(
    "Writing configuration parameters to: ",
    config_file_path,
    "\n"
  )
  write(jsonlite::toJSON(config, pretty = TRUE), config_file_path)

  return(config_file_path)
}

#' @importFrom Matrix tail
#' @importFrom stringr str_split
#' @importFrom jsonlite fromJSON
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

  # argument verificiation - use load_config for backward compatibility
  config <- load_config(config_file)

  if (config$isoform_parameters$downsample_ratio > 1 ||
        config$isoform_parameters$downsample_ratio <= 0) {
    stop("downsample_ratio should be between 0 and 1")
  }
  if (!is.null(fastq) &&
        any(!file.exists(fastq))) {
    stop(paste0("Make sure ", fastq, " exists."))
  }
  if (!file.exists(annotation)) {
    stop(paste0("Make sure ", annotation, " exists."))
  }
  if (!file.exists(genome_fa)) {
    stop(paste0("Make sure ", genome_fa, " exists."))
  }

  if (!is.null(genome_bam)) {
    if (any(!file.exists(genome_bam))) {
      stop("Make sure genome_bam exists")
    }
  }

  if (config$pipeline_parameters$bambu_isoform_identification) {
    if (!(stringr::str_ends(annotation, ".gtf") ||
            stringr::str_ends(annotation, ".gtf.gz"))) {
      stop("Bambu requires GTF format for annotation file.\n")
    }
  }

  if (config$pipeline_parameters$do_transcript_quantification &&
        config$pipeline_parameters$oarfish_quantification &&
        !config$pipeline_parameters$do_gene_quantification) {
    warning(
      "You have set to use oarfish quantification without gene quantification.",
      " Oarfish currently does not collapse UMIs, and gene quantification ",
      "performs UMI collapsing. You may want to set do_gene_quantification to ",
      "TRUE for more accurate results."
    )
  }

  return(list(config = config))
}
