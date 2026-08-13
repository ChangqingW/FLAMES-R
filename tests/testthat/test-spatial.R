test_that("create_spe errors on non-SingleCellExperiment input", {
  expect_error(create_spe(list(),       spatial_barcode_file = "dummy.txt"), "SingleCellExperiment")
  expect_error(create_spe(data.frame(), spatial_barcode_file = "dummy.txt"), "SingleCellExperiment")
  expect_error(create_spe("string",     spatial_barcode_file = "dummy.txt"), "SingleCellExperiment")
})

test_that("create_spe builds a SpatialExperiment from an SCE + barcode file + manual alignment", {
  barcodes <- c("AAAC", "AAAG", "AAAT")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0L, nrow = 2, ncol = 3,
                                  dimnames = list(c("g1", "g2"), barcodes)))
  )
  dir <- withr::local_tempdir()

  # spatial barcode file: `barcode col row`, whitespace separated, no header
  bc_file <- file.path(dir, "barcodes.txt")
  writeLines(c("AAAC 1 1", "AAAG 2 1", "AAAT 1 2"), bc_file)

  # manual alignment json supplies the imageX/imageY spatial coordinates.
  # row/col are 0-based here and get +1'd internally to match the barcode file.
  json <- file.path(dir, "align.json")
  writeLines(
    jsonlite::toJSON(
      list(oligo = data.frame(
        row = c(0L, 0L, 1L), col = c(0L, 1L, 0L),
        imageX = c(10, 30, 50), imageY = c(20, 40, 60),
        tissue = c(1L, 1L, 0L)
      )),
      auto_unbox = TRUE
    ),
    json
  )

  spe <- create_spe(sce, spatial_barcode_file = bc_file, mannual_align_json = json)

  expect_s4_class(spe, "SpatialExperiment")
  expect_equal(ncol(spe), 3L)
  expect_true("in_tissue" %in% names(SummarizedExperiment::colData(spe)))
  expect_setequal(colnames(spe), barcodes)
})
