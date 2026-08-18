# minimal bam file for tests
local_bam <- function(envir = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = envir)
  sam <- file.path(dir, "reads.sam")
  seq_a <- "ACGTACGTAAACGTACGTAC" # base 10 = A
  seq_t <- "ACGTACGTATACGTACGTAC" # base 10 = T
  qual <- paste(rep("I", 20), collapse = "")
  # `id` is the original read id; "1of1" is the <n>of<total> suffix flexiplex
  # appends (a single barcode was found in the read).
  read <- function(bc, umi, id, seq) {
    sprintf(
      "%s_%s#%s1of1\t0\tchr1\t1\t60\t20M\t*\t0\t0\t%s\t%s\tCB:Z:%s\tUB:Z:%s",
      bc, umi, id, seq, qual, bc, umi
    )
  }
  writeLines(c(
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:100",
    read("AAACCCAAGAAACGGT", "UMIAAAAAAA0", "read0", seq_a),
    read("AAACCCAAGAAACGGT", "UMIAAAAAAA1", "read1", seq_a),
    read("TTTGGGCCAAAACCCA", "UMITTTTTTT0", "read2", seq_t),
    read("TTTGGGCCAAAACCCA", "UMITTTTTTT1", "read3", seq_t)
  ), sam)
  Rsamtools::asBam(sam, sub("\\.sam$", "", sam),
    indexDestination = TRUE, overwrite = TRUE
  )
}

# A coordinate-sorted, indexed BAM with many reads on one transcript-like
# reference ("tx1"), enough to clear get_coverage()'s default min_counts = 10.
local_coverage_bam <- function(envir = parent.frame(), n = 15) {
  dir <- withr::local_tempdir(.local_envir = envir)
  sam <- file.path(dir, "cov.sam")
  seq <- paste(rep("A", 30), collapse = "")
  qual <- paste(rep("I", 30), collapse = "")
  reads <- vapply(seq_len(n), function(i) {
    pos <- ((i - 1) %% 40) + 1
    sprintf("r%d\t0\ttx1\t%d\t60\t30M\t*\t0\t0\t%s\t%s", i, pos, seq, qual)
  }, character(1))
  writeLines(c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:tx1\tLN:120", reads), sam)
  Rsamtools::asBam(sam, sub("\\.sam$", "", sam), indexDestination = TRUE, overwrite = TRUE)
}

# A single-gene annotation GRanges spanning chr1:1-100 (matches local_bam()),
# with the type/gene_id/gene_name mcols that find_variants* expect.
bam_gene_annotation <- function() {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100), strand = "+")
  S4Vectors::mcols(gr)$type <- "gene"
  S4Vectors::mcols(gr)$gene_id <- "G"
  S4Vectors::mcols(gr)$gene_name <- "G"
  gr
}

# An all-A reference matching local_bam()'s chr1, so every non-A base pileups as
# a variant.
bam_reference <- function() {
  Biostrings::DNAStringSet(c(chr1 = paste(rep("A", 100), collapse = "")))
}
