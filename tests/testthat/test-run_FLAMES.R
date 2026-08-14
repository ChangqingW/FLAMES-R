# test the full FLAMES pipeline for each pipeline

run_pipeline_case <- function(type, bambu, oarfish) {
  pipeline <- example_pipeline(type)
  pipeline@config$pipeline_parameters$bambu_isoform_identification <- bambu
  pipeline@config$pipeline_parameters$oarfish_quantification <- oarfish
  pipeline@controllers <- list() # covr might be failing due to this???
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

# the full combination will take very long
# should we reduce the number of senarios tested for certain CI branches?
run_flames_cases <- expand.grid(
  type = c("BulkPipeline", "SingleCellPipeline", "MultiSampleSCPipeline"),
  bambu = c(TRUE, FALSE),
  oarfish = c(TRUE, FALSE)
)

for (i in seq_len(nrow(run_flames_cases))) {
  local({
    .case <- run_flames_cases[i, ]
    test_that(
      sprintf("run_FLAMES completes: %s (bambu=%s, oarfish=%s)", .case$type, .case$bambu, .case$oarfish),
      {
        result <- run_pipeline_case(.case$type, .case$bambu, .case$oarfish)
        check_result(result)
      }
    )
  })
}
