test_that("barcode_output_file_identical", {
  outdir <- tempfile()
  dir.create(outdir)
  bc_allow <- local_bc_allow()

  find_barcode(
    max_bc_editdistance = 2, max_flank_editdistance = 8,
    fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
    barcodes_file = bc_allow,
    reads_out = file.path(outdir, "out.fq"),
    stats_out = file.path(outdir, "stats.tsv"),
    threads = 1, pattern = barcode_pattern(), TSO_seq = "", TSO_prime = 3, full_length_only = FALSE, strand = "+"
  )

  expect_identical(
    read.delim(test_path("bc_stat")),
    read.delim(file.path(outdir, "stats.tsv"))[, -1]
  )
  expect_identical(
    readLines(system.file('extdata', 'fastq', 'demultiplexed.fq.gz', package = 'FLAMES')),
    readLines(file.path(outdir, "out.fq"), n = 40)
  )
})
test_that("multiple fastq files as one sample", {
  outdirx <- tempfile()
  outdiry <- tempfile()
  fastq_dir <- tempfile()
  c(outdirx, outdiry, fastq_dir) |>
    sapply(dir.create)

  bc_allow <- local_bc_allow()
  # split the fastq file
  lines <- readLines(system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"))
  i <- 1
  while (length(lines) > 0) {
    writeLines(lines[1:400], file.path(fastq_dir, paste0("musc_rps24_", i, ".fastq")))
    lines <- lines[-(1:400)]
    i <- i + 1
  }

  x <- find_barcode(
    max_bc_editdistance = 2, max_flank_editdistance = 8,
    fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
    barcodes_file = bc_allow,
    reads_out = file.path(outdirx, "out.fq"),
    stats_out = file.path(outdirx, "stats.tsv"),
    threads = 1, pattern = barcode_pattern(), TSO_seq = "", TSO_prime = 3, full_length_only = FALSE
  )

  y <- find_barcode(
    max_bc_editdistance = 2, max_flank_editdistance = 8,
    fastq = fastq_dir,
    barcodes_file = bc_allow,
    reads_out = file.path(outdiry, "out.fq"),
    stats_out = file.path(outdiry, "stats.tsv"),
    threads = 1, pattern = barcode_pattern(), TSO_seq = "", TSO_prime = 3, full_length_only = FALSE
  )

  expect_identical(
    read.delim(file.path(outdirx, "stats.tsv")),
    read.delim(file.path(outdiry, "stats.tsv"))
  )
  expect_identical(
    readLines(file.path(outdirx, "out.fq")),
    readLines(file.path(outdiry, "out.fq"))
  )
})

test_that("reverse complement with strand = '-'", {
  outdir <- tempfile()
  dir.create(outdir)
  bc_allow <- local_bc_allow()

  find_barcode(
    max_bc_editdistance = 2, max_flank_editdistance = 8,
    fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
    barcodes_file = bc_allow,
    reads_out = file.path(outdir, "out.fq"),
    stats_out = file.path(outdir, "stats.tsv"),
    threads = 1, pattern = barcode_pattern(), TSO_seq = "", TSO_prime = 3, full_length_only = FALSE, strand = "-"
  )

  expect_identical(
    read.delim(test_path("bc_stat")),
    read.delim(file.path(outdir, "stats.tsv"))[, -1]
  )

  y <- readLines(
    system.file('extdata', 'fastq', 'demultiplexed.fq.gz', package = 'FLAMES')
  )[seq(2, 40, by = 4)] |>
    Biostrings::DNAStringSet() |>
    Biostrings::reverseComplement() |>
    as.character()

  expect_identical(
    readLines(file.path(outdir, "out.fq"), n = 40)[seq(2, 40, by = 4)],
    y
  )
})

