#' BLAZE Assign reads to cell barcodes. 
#'
#' @description
#' Uses BLAZE to generate barcode list and assign reads to cell barcodes. 
#' @param expect_cells Integer, expected number of cells. Note: this could be just a rough estimate. E.g., the targeted number of cells.
#' @param fq_in File path to the fastq file used as a query sequence file
#' @param additional_args Additional command line style arguments to be passed to BLAZE. E.g. c("--10x-kit-version", "3v3")
#' @param ... Additional BLAZE configuration parameters. E.g., setting
#'  `'output-prefix'='some_prefix'` is equivalent to specifying `--output-prefix some_prefix` in BLAZE; Similarly,
#'  `overwrite=TRUE` is equivalent to switch on the `--overwrite` option. Note that the specified parameters will
#'  override the parameters specified in the configuration file. All available options can be found at https://github.com/shimlab/BLAZE.
#'
#' @return A \code{data.frame} summarising the reads aligned. Other outputs are written to disk. 
#' The details of the output files can be found at https://github.com/shimlab/BLAZE.
#'
#' @export
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' fastq <- system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES")
#' blaze(
#'   expect_cells = 10, fastq,
#'   "output-prefix" = file.path(outdir, ""),
#'   "output-fastq" = file.path(outdir, "output.fastq"),
#'   overwrite=TRUE
#' )
#'
#' @importFrom reticulate import_from_path dict
#' @importFrom basilisk basiliskRun
#' @export
blaze <- function(expect_cells, fq_in, additional_args = NULL, ...) {
        # prepare command-line-style arguments for blaze        
        blaze_argv <- paste('--expect-cells ', expect_cells)
        
        # include additional blaze parameters
        blaze_config <- list(...)
        
        
        if (length(blaze_config) > 0) {
            # handle the switch options first as they do not have values
            if ('overwrite' %in% names(blaze_config) && blaze_config$`overwrite` == TRUE) {
                blaze_argv <- paste(blaze_argv, '--overwrite --minimal_stdout ')
            }
            
            print(blaze_config)
            if ('high-sensitivity-mode' %in% names(blaze_config)
                && blaze_config$`high-sensitivity-mode` == TRUE) {
                blaze_argv <- paste(blaze_argv, '--high-sensitivity-mode ')
            }

            # remove the switch options from the list and add the rest
            blaze_config['overwrite'] <- NULL
            blaze_config['high-sensitivity-mode'] <- NULL
            for (arg in names(blaze_config)) {
                blaze_argv <- paste(blaze_argv, paste0('--',arg), blaze_config[arg])}
        }

        blaze_argv <- paste(blaze_argv, paste0(additional_args, collapse = " "), fq_in)


        # run blaze
        ret <-
            basiliskRun(env = flames_env, fun = function(blaze_argv) {
                # blaze_path <- system.file("blaze", package = "FLAMES")
                cat("Running BLAZE...\n")
                cat("Argument: ", blaze_argv, "\n")
                blaze <-
                    reticulate::import("blaze")
                ret <-
                    blaze$blaze(blaze_argv)

                ret
            }, blaze_argv = blaze_argv
            )
        #ret # return the filename of demultiplexed fastq
    }
