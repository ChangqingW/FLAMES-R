test_that("sc_mutations errors when seqnames and positions differ in length", {
  expect_error(
    sc_mutations("dummy.bam", seqnames = c("chr1", "chr2"), positions = c(100)),
    "seqnames not the same length as positions"
  )
})

test_that("create_spe errors immediately when sce is not a SingleCellExperiment", {
  expect_error(create_spe(list(),      spatial_barcode_file = "dummy.txt"), "SingleCellExperiment")
  expect_error(create_spe(data.frame(), spatial_barcode_file = "dummy.txt"), "SingleCellExperiment")
  expect_error(create_spe("string",    spatial_barcode_file = "dummy.txt"), "SingleCellExperiment")
})

