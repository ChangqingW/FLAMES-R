make_pipeline <- function(type, bambu, oarfish, controller, demux = "flexiplex") {
  pipeline <- example_pipeline(type)
  pipeline@config$pipeline_parameters$bambu_isoform_identification <- bambu
  pipeline@config$pipeline_parameters$oarfish_quantification <- oarfish
  pipeline@config$pipeline_parameters$demultiplexer <- demux
  if (demux == "BLAZE" && methods::is(pipeline, "FLAMES.SingleCellPipeline")) {
    pipeline@expect_cell_number <- rep(100L, length(pipeline@fastq))
  }
  pipeline@controllers <- if (controller) {
    list(default = crew::crew_controller_local())
  } else {
    list()
  }
  pipeline
}

run_pipeline_case <- function(type, bambu, oarfish, controller) {
  run_FLAMES(make_pipeline(type, bambu, oarfish, controller))
}

check_result <- function(result) {
  # elavate failed steps
  if (length(result@last_error) > 0) {
    fail(sprintf(
      "pipeline step '%s' failed: %s",
      result@last_error$step, conditionMessage(result@last_error$error)
    ))
    return(invisible())
  }
  expect_true(all(result@completed_steps[result@steps]))
  exp <- experiment(result)
  expect_false(is.null(exp))
  if (is.list(exp) && !methods::is(exp, "SummarizedExperiment")) {
    # multi-sample pipelines return a list of SingleCellExperiment objects
    expect_gt(length(exp), 0)
    expect_gt(nrow(exp[[1]]), 0)
  } else {
    expect_gt(nrow(exp), 0)
    expect_gt(ncol(exp), 0)
  }
}

run_flames_cases <- data.frame(
  type = c(
    "BulkPipeline", "SingleCellPipeline", "MultiSampleSCPipeline",
    "BulkPipeline", "SingleCellPipeline", "MultiSampleSCPipeline"
  ),
  bambu      = c(TRUE,  FALSE, FALSE, FALSE, FALSE, FALSE),
  oarfish    = c(TRUE,  FALSE, TRUE,  FALSE, FALSE, FALSE),
  controller = c(FALSE, FALSE, FALSE, TRUE,  TRUE,  TRUE),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(run_flames_cases))) {
  local({
    .case <- run_flames_cases[i, ]
    test_that(
      sprintf(
        "run_FLAMES completes: %s (bambu=%s, oarfish=%s, controller=%s)",
        .case$type, .case$bambu, .case$oarfish, .case$controller
      ),
      {
        result <- run_pipeline_case(.case$type, .case$bambu, .case$oarfish, .case$controller)
        check_result(result)
      }
    )
  })
}

# cover BLAZE as well. just the demultiplex step.
# BLAZE inside local controller will crash covr
for (i in c("SingleCellPipeline", "MultiSampleSCPipeline")) {
  local({
    test_that(
      sprintf("barcode_demultiplex completes with BLAZE: %s (no controller)", i),
      {
        pipeline <- make_pipeline(i, bambu = FALSE, oarfish = FALSE,
                                  controller = FALSE, demux = "BLAZE")
        pipeline <- run_step(pipeline, "barcode_demultiplex", disable_controller = FALSE)
        expect_true(pipeline@completed_steps[["barcode_demultiplex"]])
        expect_true(all(nzchar(pipeline@demultiplexed_fastq)))
        expect_true(all(file.exists(pipeline@demultiplexed_fastq)))
      }
    )
  })
}
