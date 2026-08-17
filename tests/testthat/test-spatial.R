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

# A minimal SpatialExperiment carrying a 2x2 fake tissue image + imageX/imageY
# coords + a counts assay, enough to exercise the plotting engine (magick etc.).
make_image_spe <- function() {
  counts <- matrix(c(5, 1, 0, 2, 3, 4), nrow = 2, byrow = TRUE,
                   dimnames = list(c("f1", "f2"), paste0("spot", 1:3)))
  coords <- matrix(c(1, 1, 2, 1, 1, 2), ncol = 2, byrow = TRUE,
                   dimnames = list(NULL, c("imageX", "imageY")))
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts), spatialCoords = coords
  )
  ras <- grDevices::as.raster(matrix(c("#FF0000", "#00FF00", "#0000FF", "#FFFFFF"), 2, 2))
  spi <- SpatialExperiment::SpatialImage(ras)
  SpatialExperiment::imgData(spe) <- S4Vectors::DataFrame(
    sample_id = "sample01", image_id = "hires", data = I(list(spi)), scaleFactor = 1
  )
  spe
}

test_that("plot_spatial_feature returns a ggplot over the tissue image", {
  spe <- make_image_spe()
  expect_s3_class(plot_spatial_feature(spe, "f1"), "ggplot")
  # grayscale = FALSE path + numeric-vector feature of length ncol(spe)
  expect_s3_class(
    plot_spatial_feature(spe, feature = c(1, 2, 3), grayscale = FALSE),
    "ggplot"
  )
})

test_that("plot_spatial_feature rejects a wrongly-sized feature", {
  spe <- make_image_spe()
  expect_error(plot_spatial_feature(spe, feature = c(1, 2)), "length 1 or ncol")
})

test_that("plot_spatial_pie returns a ggplot for multiple features", {
  spe <- make_image_spe()
  expect_s3_class(plot_spatial_pie(spe, c("f1", "f2")), "ggplot")
})

test_that("plot_spatial errors when the SpatialExperiment has no image data", {
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = matrix(1, 1, 3, dimnames = list("f1", paste0("s", 1:3)))),
    spatialCoords = matrix(c(1, 1, 2, 1, 1, 2), ncol = 2, byrow = TRUE,
                           dimnames = list(NULL, c("imageX", "imageY")))
  )
  expect_error(FLAMES:::plot_spatial(spe, opacity = 50, feature = 1:3), "No image data")
})
