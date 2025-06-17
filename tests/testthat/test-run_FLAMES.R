test_that("run_FLAMES() completes and returns experiment for all pipeline types and option combinations", {
  pipeline_types <- c("BulkPipeline", "SingleCellPipeline", "MultiSampleSCPipeline")
  bools <- c(TRUE, FALSE)
  controllers <- list(
    list(),
    list(
      genome_alignment = crew::crew_controller_local(),
      read_realignment = crew::crew_controller_local()
    )
  )

  for (ptype in pipeline_types) {
    for (bambu in bools) {
      for (oarfish in bools) {
        for (controller in controllers) {
          test_name <- sprintf(
            "run_FLAMES(%s, bambu=%s, oarfish=%s, controller=%s)",
            ptype, bambu, oarfish,
            ifelse(inherits(controller, "crew_class_controller"), "crew_controller_local", "none")
          )
          message("Testing: ", test_name)

          pipeline <- example_pipeline(ptype)
          pipeline@config$pipeline_parameters$bambu_isoform_identification <- bambu
          pipeline@config$pipeline_parameters$oarfish_quantification <- oarfish
          pipeline@controllers <- controller

          result <- run_FLAMES(pipeline)
          expect_false(
            is.null(experiment(result)),
            info = paste("experiment is NULL for", test_name)
          )
        }
      }
    }
  }
})
