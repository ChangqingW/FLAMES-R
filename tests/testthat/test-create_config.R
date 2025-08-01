test_that("create_config check (legacy test)", {
  outdir <- tempfile()
  dir.create(outdir)
  # Test that the refactored create_config still handles old-style parameters
  test_conf <- create_config(outdir, max_bc_editdistance = 123)
  expect_equal(jsonlite::fromJSON(test_conf)$barcode_parameters$max_bc_editdistance, 123)

  test_lists <- jsonlite::fromJSON(system.file("extdata/config_sclr_nanopore_3end.json", package = "FLAMES"))
  test_lists$barcode_parameters$max_bc_editdistance <- as.integer(123)
  test_lists$comment <- NULL
  expect_identical(test_lists, jsonlite::fromJSON(test_conf))
})
