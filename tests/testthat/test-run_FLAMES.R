run_pipeline_case <- function(type, bambu, oarfish, controller) {
  pipeline <- example_pipeline(type)
  pipeline@config$pipeline_parameters$bambu_isoform_identification <- bambu
  pipeline@config$pipeline_parameters$oarfish_quantification <- oarfish
  pipeline@controllers <- if (controller) {
    list(default = crew::crew_controller_local())
  } else {
    list()
  }
  run_FLAMES(pipeline)
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
