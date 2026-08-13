# no external aligners or python are involved

test_that("get_GRangesList imports a GTF path into a per-transcript GRangesList", {
  gtf <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
  res <- FLAMES:::get_GRangesList(gtf)

  expect_named(res, c("grl", "rowdata"))
  expect_s4_class(res$grl, "CompressedGRangesList")
  expect_equal(length(res$grl), 10L)          # 10 transcripts in Rps24
  expect_equal(nrow(res$rowdata), 10L)
  expect_true("transcript_id" %in% colnames(res$rowdata))
  # rowdata is aligned row-for-row with the GRangesList
  expect_identical(names(res$grl), res$rowdata$transcript_id)
})

test_that("get_GRangesList accepts a GRanges and splits it by transcript_id", {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1, 100, 1), end = c(50, 150, 50)),
    strand = "+"
  )
  S4Vectors::mcols(gr)$transcript_id <- c("txA", "txA", "txB")
  S4Vectors::mcols(gr)$type <- "exon"

  res <- FLAMES:::get_GRangesList(gr)
  expect_equal(length(res$grl), 2L)
  expect_setequal(names(res$grl), c("txA", "txB"))
  expect_equal(nrow(res$rowdata), 2L)
})

test_that("get_GRangesList returns a GRangesList input unchanged with NULL rowdata", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10))
  grl <- GenomicRanges::GRangesList(txA = gr)
  res <- FLAMES:::get_GRangesList(grl)
  expect_identical(res$grl, grl)
  expect_null(res$rowdata)
})

test_that("fake_stranded_gff returns the input unchanged when all features are stranded", {
  gtf <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
  expect_identical(FLAMES:::fake_stranded_gff(gtf), gtf)
})

test_that("fake_stranded_gff rewrites '.' strands to '+' and warns", {
  dir <- withr::local_tempdir()
  gtf <- file.path(dir, "unstranded.gtf")
  writeLines(c(
    'chr1\ttest\texon\t1\t100\t.\t.\t.\tgene_id "g1"; transcript_id "t1";',
    'chr1\ttest\texon\t200\t300\t.\t.\t.\tgene_id "g1"; transcript_id "t1";'
  ), gtf)

  expect_warning(out <- FLAMES:::fake_stranded_gff(gtf), "not stranded")
  expect_false(identical(out, gtf))
  expect_true(file.exists(out))
  d <- read.csv(out, sep = "\t", header = FALSE, comment.char = "#")
  expect_true(all(d[, 7] == "+"))
})

test_that("annotation_to_fasta writes a transcriptome FASTA and its .fai index", {
  gtf <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
  fa <- system.file("extdata", "rps24.fa.gz", package = "FLAMES")
  out <- file.path(withr::local_tempdir(), "transcriptome.fa")

  annotation_to_fasta(gtf, fa, out)

  expect_true(file.exists(out))
  expect_true(file.exists(paste0(out, ".fai")))
  seqs <- Biostrings::readDNAStringSet(out)
  expect_gt(length(seqs), 0)
  expect_true(all(Biostrings::width(seqs) > 0))
})
