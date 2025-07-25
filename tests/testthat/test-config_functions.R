test_that("merge_configs_recursive works correctly", {
  # Test basic merging
  default <- list(
    a = 1,
    b = list(x = 10, y = 20),
    c = list(nested = list(deep = "default"))
  )
  
  user <- list(
    a = 2,  # override
    b = list(y = 30),  # partial override
    d = "new"  # new parameter
  )
  
  result <- merge_configs_recursive(default, user)
  
  expect_equal(result$a, 2)  # user override
  expect_equal(result$b$x, 10)  # default preserved
  expect_equal(result$b$y, 30)  # user override
  expect_equal(result$c$nested$deep, "default")  # default preserved
  expect_equal(result$d, "new")  # new parameter added
  
  # Test with NULL user config
  result_null <- merge_configs_recursive(default, NULL)
  expect_equal(result_null, default)
  
  # Test comment removal
  default_with_comment <- list(comment = "test", a = 1)
  user_simple <- list(a = 2)
  result_no_comment <- merge_configs_recursive(default_with_comment, user_simple)
  expect_false("comment" %in% names(result_no_comment))
  expect_equal(result_no_comment$a, 2)
})

test_that("create_config with dot notation works", {
  outdir <- tempfile()
  dir.create(outdir)
  
  # Test basic functionality (backward compatibility)
  test_conf <- create_config(outdir, max_bc_editdistance = 123)
  config <- jsonlite::fromJSON(test_conf)
  expect_equal(config$barcode_parameters$max_bc_editdistance, 123)
  
  # Test dot notation for nested parameters
  test_conf_nested <- create_config(outdir, 
    threads = 16,
    "barcode_parameters.max_bc_editdistance" = 5,
    "barcode_parameters.pattern.primer" = "TESTSEQ",
    "isoform_parameters.min_sup_cnt" = 10
  )
  
  config_nested <- jsonlite::fromJSON(test_conf_nested)
  expect_equal(config_nested$pipeline_parameters$threads, 16)
  expect_equal(config_nested$barcode_parameters$max_bc_editdistance, 5)
  expect_equal(config_nested$barcode_parameters$pattern$primer, "TESTSEQ")
  expect_equal(config_nested$isoform_parameters$min_sup_cnt, 10)
  
  # Test SIRV type
  test_conf_sirv <- create_config(outdir, type = "SIRV")
  config_sirv <- jsonlite::fromJSON(test_conf_sirv)
  expect_true(config_sirv$alignment_parameters$no_flank)
  
  # Test error on unnamed parameters in ... 
  outdir2 <- tempfile()
  dir.create(outdir2)
  expect_error({
    create_config(outdir2, "sc_3end", 123)
  }, "Parameters must be named")
  
  # Test error on invalid type
  expect_error(create_config(outdir, type = "invalid"), "Unrecognised config type")
})

test_that("load_config provides backward compatibility", {
  outdir <- tempfile()
  dir.create(outdir)
  
  # Create a minimal config file (simulating old version)
  minimal_config <- list(
    pipeline_parameters = list(threads = 4),
    barcode_parameters = list(max_bc_editdistance = 3)
  )
  
  config_file <- file.path(outdir, "minimal_config.json")
  write(jsonlite::toJSON(minimal_config, pretty = TRUE), config_file)
  
  # Load with load_config - should fill in missing defaults
  loaded_config <- load_config(config_file)
  
  # Check that user values are preserved
  expect_equal(loaded_config$pipeline_parameters$threads, 4)
  expect_equal(loaded_config$barcode_parameters$max_bc_editdistance, 3)
  
  # Check that missing defaults are filled in
  expect_true("do_genome_alignment" %in% names(loaded_config$pipeline_parameters))
  expect_true("pattern" %in% names(loaded_config$barcode_parameters))
  expect_true("isoform_parameters" %in% names(loaded_config))
  
  # Check specific default values are present
  expect_equal(loaded_config$pipeline_parameters$do_genome_alignment, TRUE)
  expect_equal(loaded_config$isoform_parameters$max_dist, 10)
})

test_that("set_nested_param works correctly", {
  config <- list(
    a = 1,
    b = list(x = 10)
  )
  
  # Test setting existing nested parameter
  result1 <- set_nested_param(config, "b.x", 20)
  expect_equal(result1$b$x, 20)
  
  # Test creating new nested structure
  result2 <- set_nested_param(config, "c.d.e", "new_value")
  expect_equal(result2$c$d$e, "new_value")
  
  # Test setting new parameter in existing structure
  result3 <- set_nested_param(config, "b.y", 30)
  expect_equal(result3$b$y, 30)
  expect_equal(result3$b$x, 10)  # original preserved
})