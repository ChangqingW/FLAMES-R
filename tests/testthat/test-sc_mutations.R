test_that("sc_mutations errors when seqnames and positions differ in length", {
  expect_error(
    sc_mutations("dummy.bam", seqnames = c("chr1", "chr2"), positions = c(100)),
    "seqnames not the same length as positions"
  )
})

test_that("extract_nt extracts single nucleotides by seqname and position", {
  ref <- Biostrings::DNAStringSet(c(chrA = "ACGTACGT", chrB = "TTTTGGGG"))
  res <- FLAMES:::extract_nt(ref, c("chrA", "chrA", "chrB"), c(1, 3, 5))
  expect_equal(unname(res), c("A", "G", "G"))
})

test_that("homopolymer_pct is 1 within a homopolymer run and low in mixed sequence", {
  ref <- Biostrings::DNAStringSet(c(
    poly = "TAAAAAAAAAAT",   # long A run
    mix  = "ACGTACGTACGT"
  ))
  # interior position inside the A-run, flanking bases are all A
  expect_equal(unname(FLAMES:::homopolymer_pct(ref, "poly", 6, n = 3)), 1)
  # mixed neighbourhood
  expect_lt(unname(FLAMES:::homopolymer_pct(ref, "mix", 6, n = 3)), 0.5)
})

test_that("edgecase: homopolymer_pct treats terminal positions as homopolymers (returns 1)", {
  ref <- Biostrings::DNAStringSet(c(mix = "ACGTACGTACGT"))
  res <- FLAMES:::homopolymer_pct(ref, c("mix", "mix"), c(1, 12), n = 3)
  expect_equal(as.numeric(res), c(1, 1))
})

# A single gene on chr1 spanning 1-100 with two exons (1-20, 51-80);
# total exonic length is 20 + 30 = 50.
make_gene_annotation <- function(strand = "+", gene_id = "G1") {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1, 1, 51), end = c(100, 20, 80)),
    strand = strand
  )
  S4Vectors::mcols(gr)$type <- c("gene", "exon", "exon")
  S4Vectors::mcols(gr)$gene_id <- gene_id
  S4Vectors::mcols(gr)$gene_name <- "TestGene"
  gr
}

test_that("mutation_positions_single computes relative positions, NA for intronic sites", {
  annot <- make_gene_annotation("+")
  mut <- data.frame(seqnames = "chr1", pos = c(10, 60, 30), region = "TestGene")
  res <- FLAMES:::mutation_positions_single(mut, annot, type = "relative", verbose = FALSE)
  expect_equal(res, c(0.2, 0.6, NA))
})

test_that("mutation_positions_single supports TSS and TES distances", {
  annot <- make_gene_annotation("+")
  mut <- data.frame(seqnames = "chr1", pos = c(10, 60), region = "TestGene")
  expect_equal(FLAMES:::mutation_positions_single(mut, annot, type = "TSS", verbose = FALSE), c(10, 30))
  expect_equal(FLAMES:::mutation_positions_single(mut, annot, type = "TES", verbose = FALSE), c(40, 20))
})

test_that("mutation_positions_single flips coordinates on the negative strand", {
  annot <- make_gene_annotation("-")
  mut <- data.frame(seqnames = "chr1", pos = c(10, 60), region = "TestGene")
  expect_equal(
    FLAMES:::mutation_positions_single(mut, annot, type = "relative", verbose = FALSE),
    c(0.8, 0.4)
  )
})

test_that("mutation_positions_single errors on multiple regions or multiple genes", {
  annot <- make_gene_annotation("+")
  mut_multi <- data.frame(seqnames = "chr1", pos = c(10, 60), region = c("A", "B"))
  expect_error(
    FLAMES:::mutation_positions_single(mut_multi, annot, type = "relative", verbose = FALSE),
    "region"
  )
  annot2 <- make_gene_annotation("+", gene_id = c("G1", "G1", "G2"))
  mut <- data.frame(seqnames = "chr1", pos = 10, region = "TestGene")
  expect_error(
    FLAMES:::mutation_positions_single(mut, annot2, type = "relative", verbose = FALSE),
    "gene_id"
  )
})

test_that("mutation_positions_single errors when a mutation is outside the gene", {
  annot <- make_gene_annotation("+")
  mut <- data.frame(seqnames = "chr1", pos = c(10, 200), region = "TestGene")
  expect_error(
    FLAMES:::mutation_positions_single(mut, annot, type = "relative", verbose = FALSE),
    "not within the gene region"
  )
})

test_that("mutation_positions returns relative positions for the bundled Rps24 annotation", {
  variants <- data.frame(
    seqnames = rep("chr14", 8),
    pos = c(1084, 1085, 1217, 1384, 2724, 2789, 5083, 5147),
    region = rep("Rps24", 8)
  )
  gtf <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
  res <- mutation_positions(mutations = variants, annotation = gtf)
  expect_length(res, 8)
  expect_true(all(is.na(res) | (res >= 0 & res <= 1)))
})

test_that("mutation_positions with bin = TRUE returns a 100-column matrix per gene", {
  variants <- data.frame(
    seqnames = rep("chr14", 8),
    pos = c(1084, 1085, 1217, 1384, 2724, 2789, 5083, 5147),
    region = rep("Rps24", 8)
  )
  gtf <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
  res <- mutation_positions(mutations = variants, annotation = gtf, bin = TRUE)
  expect_true(is.matrix(res))
  expect_equal(ncol(res), 100L)
  expect_equal(rownames(res), "Rps24")
})

test_that("mutation_positions errors on an unknown type", {
  variants <- data.frame(seqnames = "chr14", pos = 1084, region = "Rps24")
  gtf <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
  expect_error(
    mutation_positions(variants, gtf, type = "bad"),
    "'type' must be one of"
  )
})

test_that("sc_genotype classifies barcodes as alt / ref / NA by thresholds", {
  snps_tb <- tibble::tibble(
    seqname = "chr1",
    pos = 100L,
    barcode = c("bc1", "bc1", "bc2", "bc3"),
    allele = c("A", "G", "G", "A"),
    allele_count = c(5, 1, 4, 1),
    pct = c(5 / 6, 1 / 6, 1.0, 0.05),
    cell_total_reads = c(6, 6, 4, 20)
  )
  gt <- sc_genotype(
    snps_tb, ref = "G", alt = "A", seqname = "chr1", pos = 100,
    alt_min_count = 2, alt_min_pct = 0.1, ref_min_count = 1, ref_min_pct = 1
  )

  expect_s3_class(gt, "tbl_df")
  expect_true(all(c("barcode", "allele_count_ref", "pct_ref",
                    "allele_count_alt", "pct_alt", "genotype") %in% names(gt)))

  g <- setNames(gt$genotype, gt$barcode)
  expect_equal(g[["bc1"]], "alt")
  expect_equal(g[["bc2"]], "ref")
  expect_true(is.na(g[["bc3"]]))
})

test_that("sc_plot_genotype returns a ggplot on a reduced-dim embedding", {
  barcodes <- paste0("bc", 1:6)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0, nrow = 2, ncol = 6,
                                  dimnames = list(c("g1", "g2"), barcodes)))
  )
  umap <- matrix(as.numeric(1:12), ncol = 2,
                 dimnames = list(barcodes, c("UMAP1", "UMAP2")))
  SingleCellExperiment::reducedDims(sce)$UMAP <- umap

  genotype_tb <- tibble::tibble(
    barcode = c("bc1", "bc2", "bc3"),
    genotype = c("alt", "ref", "alt")
  )
  p <- sc_plot_genotype(sce, genotype_tb)
  expect_s3_class(p, "ggplot")
})
