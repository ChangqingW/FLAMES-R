#' @importFrom tibble as_tibble tibble
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect matches
#' @importFrom dplyr group_by mutate ungroup
variant_count_tb <- function(bam_path, seqname, pos, indel, verbose = TRUE) {
  # allele by barcode matrix (value: read count)
  tryCatch(
    {
      tb <- variant_count_matrix( # throws Rcpp::exception when no reads at pos
        bam_path = bam_path,
        seqname = seqname, pos = pos, indel = indel, verbose = verbose
      ) |>
        tibble::as_tibble(rownames = "allele")
      if (ncol(tb) <= 1) { # no cols besides "allele"
        message(sprintf("No reads found at %s:%d in %s (All RefSkips)", seqname, pos, bam_path))
        return(tibble::tibble())
      }
        # pivot to long format: allele, barcode, allele_count
      tb |>
        tidyr::pivot_longer(
          cols = -tidyselect::matches("allele"),
          values_to = "allele_count", names_to = "barcode"
        ) |>
        dplyr::group_by(barcode) |>
        dplyr::mutate(
          cell_total_reads = sum(allele_count),
          pct = allele_count / cell_total_reads,
          pos = pos, seqname = seqname
        ) |>
        dplyr::ungroup()
    },
    error = function(e) {
      if (inherits(e, "Rcpp::exception") & conditionMessage(e) == "Failed to fetch an alignment") {
        message(paste0("No reads found at ", seqname, ":", pos, " in ", bam_path))
        message("Returning empty tibble")
        return(tibble::tibble())
      } else {
        stop(e)
      }
    }
  )
}

#' Variant count for single-cell data
#'
#' Count the number of reads supporting each variants at the given positions for each cell.
#'
#' @importFrom BiocParallel bplapply bpmapply MulticoreParam
#' @importFrom dplyr mutate select bind_rows
#'
#' @param bam_path character(1) or character(n): path to the bam file(s) aligned to the
#' reference genome (NOT the transcriptome! Unless the postions are also from the transcriptome).
#' @param seqnames character(n): chromosome names of the postions to count alleles.
#' @param positions integer(n): positions, 1-based, same length as seqnames. The positions to count alleles.
#' @param indel logical(1): whether to count indels (TRUE) or SNPs (FALSE).
#' @param threads integer(1): number of threads to use. Maximum number of threads is
#' the number of bam files * number of positions.
#' @return A tibble with columns: allele, barcode, allele_count, cell_total_reads, pct, pos, seqname.
#' @examples
#' ppl <- example_pipeline("SingleCellPipeline")
#' ppl <- run_step(ppl, "barcode_demultiplex")
#' ppl <- run_step(ppl, "genome_alignment")
#' snps_tb <- sc_mutations(
#'   bam_path = ppl@genome_bam,
#'   seqnames = c("chr14", "chr14"),
#'   positions = c(1260, 2714), # positions of interest
#'   indel = FALSE
#' )
#' head(snps_tb)
#' snps_tb |>
#'   dplyr::filter(pos == 1260) |>
#'   dplyr::group_by(allele) |>
#'   dplyr::summarise(count = sum(allele_count)) # should be identical to samtools pileup
#' @export
sc_mutations <- function(bam_path, seqnames, positions, indel = FALSE, threads = 1) {
  stopifnot(
    "seqnames not the same length as positions" =
      length(seqnames) == length(positions)
  )

  if (length(bam_path) == 1) {
    # single bam file, parallelize over positions
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Got 1 bam file, parallelizing over each position ..."))
    variants <- tryCatch({
      BiocParallel::bpmapply(
        function(seqname, pos) {
          variant_count_tb(bam_path, seqname, pos, indel, verbose = FALSE)
        },
        seqname = seqnames, pos = positions, SIMPLIFY = FALSE,
        BPPARAM = BiocParallel::MulticoreParam(
          workers = threads, stop.on.error = TRUE, progressbar = TRUE
        )
      )
    }, error = identity)
  } else {
    # multiple bam files, parallelize over bam files
    # data frame of all combinations between (seqname, pos) and (bam_path, barcodes)
    args_grid <- expand.grid(
      mutation_index = seq_along(positions),
      bam_index = seq_along(bam_path),
      stringsAsFactors = FALSE
    ) |>
      dplyr::mutate(
        seqname = seqnames[mutation_index],
        pos = positions[mutation_index],
        sample_bam = bam_path[bam_index]
      ) |>
      dplyr::select(-mutation_index, -bam_index)

    message(paste0(format(Sys.time(), "%H:%M:%S "), "Multi-threading over bam files x positions ..."))
    variants <- tryCatch({
      BiocParallel::bpmapply(
        function(sample_bam, seqname, pos) {
          variant_count_tb(sample_bam, seqname, pos, indel, verbose = FALSE) |>
            dplyr::mutate(bam_file = sample_bam)
        },
        sample_bam = args_grid$sample_bam, seqname = args_grid$seqname,
        pos = args_grid$pos, SIMPLIFY = FALSE,
        BPPARAM = BiocParallel::MulticoreParam(
          workers = threads, stop.on.error = TRUE, progressbar = TRUE
        )
      )
    }, error = identity)
  }

  if (inherits(variants, "error")) {
    warning("Error occurred in `sc_mutations`, returning error object")
    return(variants)
  }

  message(paste0(format(Sys.time(), "%H:%M:%S "), "Merging results ..."))
  variants <- variants |>
    dplyr::bind_rows()
  return(variants)
}


extract_nt <- function(ref, seqname, pos) {
  mapply(function(seqname, pos) {
    as.character(ref[[seqname]][pos])
  }, seqname, pos)
}

#' @importFrom BiocParallel bplapply bpmapply MulticoreParam
#' @importFrom dplyr mutate pull
homopolymer_pct <- function(ref, seqname, pos, include_alt = FALSE, n = 3, threads = 1) {
  pcts <- tryCatch({
    BiocParallel::bpmapply(
      function(seqname, pos, include_alt) {
        if (pos == 1 | pos == length(ref[[seqname]])) {
          return(TRUE) # variant at the ends should not be considered
        }
        start <- max(1, pos - n)
        end <- min(length(ref[[seqname]]), pos + n)
        if (include_alt) {
          ref[[seqname]][start:end] |>
            as.character() |>
            strsplit("") |>
            unlist() |>
            table() |>
            as.data.frame() |>
            dplyr::mutate(pct = Freq / sum(Freq)) |>
            dplyr::pull(pct) |>
            max()
        } else {
          # exclude the position itself
          ref[[seqname]][c(start:(pos - 1), (pos + 1):end)] |>
            as.character() |>
            strsplit("") |>
            unlist() |>
            table() |>
            as.data.frame() |>
            dplyr::mutate(pct = Freq / sum(Freq)) |>
            dplyr::pull(pct) |>
            max()
        }
      },
      seqname, pos, include_alt, SIMPLIFY = TRUE,
      BPPARAM = BiocParallel::MulticoreParam(
        workers = threads, stop.on.error = TRUE, progressbar = TRUE
      )
    )
  }, error = identity)
  if (inherits(pcts, "error")) {
    warning("Error occurred while running `homopolymer_pct`, returning error object")
  }
  return(pcts)
}

# WIP: too much sequencing errrors / splice sites
# find variants in a single grange
#' @importFrom Rsamtools pileup PileupParam ScanBamParam
#' @importFrom dplyr select mutate group_by ungroup filter
#' @importFrom S4Vectors mcols
find_variants_grange <- function(bam_path, reference, gene_grange, min_nucleotide_depth,
                                 names_from) {
  # read bam file
  mutations <- Rsamtools::pileup(bam_path,
    pileupParam = Rsamtools::PileupParam(
      max_depth = .Machine$integer.max - 1, min_base_quality = 0, min_mapq = 0,
      min_nucleotide_depth = min_nucleotide_depth, min_minor_allele_depth = 0,
      distinguish_strands = FALSE, distinguish_nucleotides = TRUE,
      ignore_query_Ns = TRUE, include_deletions = TRUE, include_insertions = TRUE,
      left_bins = NULL, query_bins = NULL, cycle_bins = NULL
    ),
    scanBamParam = Rsamtools::ScanBamParam(which = gene_grange)
  ) |>
    dplyr::select(-which_label) |>
    # move insertion postions to the base at which insertion start
    # like in IGV
    # don't count insertion in sum
    dplyr::mutate(
      pos = ifelse(nucleotide == "+", pos + 1, pos),
      counts_no_ins = ifelse(nucleotide != "+", count, 0)
    ) |>
    dplyr::group_by(seqnames, pos) |>
    dplyr::mutate(sum = sum(counts_no_ins)) |>
    dplyr::ungroup() |>
    dplyr::select(-counts_no_ins) |>
    dplyr::mutate(
      freq = count / sum,
      ref = factor(extract_nt(reference, seqnames, pos))
    ) |>
    dplyr::filter(as.character(nucleotide) != as.character(ref))

  if (nrow(mutations) == 0) {
    return(mutations)
  } else {
    mutations$bam_path <- bam_path
    if (names_from %in% colnames(S4Vectors::mcols(gene_grange))) {
      mutations$region <- S4Vectors::mcols(gene_grange)[, names_from]
    } else {
      mutations$region <- NA # no gene name / gap
    }
    return(mutations)
  }
}

#' bulk variant identification
#'
#' Treat each bam file as a bulk sample and identify variants against the reference
#'
#' Each bam file is treated as a bulk sample to perform pileup and identify variants.
#' You can run \code{sc_mutations} with the variants identified with this function
#' to get single-cell allele counts. Note that reference genome FASTA files may have
#' the chromosome names field as `>chr1 1` instead of `>chr1`. You may need to remove
#' the trailing number to match the chromosome names in the bam file, for example with
#' \code{names(ref) <- sapply(names(ref), function(x) strsplit(x, " ")[[1]][1])}.
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom rtracklayer import
#' @importFrom S4Vectors mcols mcols<- DataFrame
#' @importFrom GenomeInfoDb seqnames seqlengths seqinfo
#' @importFrom GenomicRanges disjoin gaps
#' @importFrom BiocParallel bplapply bpmapply MulticoreParam
#' @importFrom dplyr bind_rows mutate
#'
#' @param bam_path character(1) or character(n): path to the bam file(s) aligned to the
#' reference genome (NOT the transcriptome!).
#' @param reference DNAStringSet: the reference genome
#' @param annotation GRanges: the annotation of the reference genome. You can load
#' a GTF/GFF annotation file with \code{anno <- rtracklayer::import(file)}.
#' @param min_nucleotide_depth integer(1): minimum read depth for a position to be
#' considered a variant.
#' @param threads integer(1): number of threads to use. Threading is done over each
#' annotated region and (if \code{annotated_region_only = FALSE}) unannotated gaps for
#' each bam file.
#' @param homopolymer_window integer(1): the window size to calculate the homopolymer
#' percentage. The homopolymer percentage is calculated as the percentage of the most
#' frequent nucleotide in a window of \code{-homopolymer_window} to \code{homopolymer_window}
#' nucleotides around the variant position, excluding the variant position itself.
#' Calculation of the homopolymer percentage is skipped when \code{homopolymer_window = 0}.
#' This is useful for filtering out Nanopore sequencing errors in homopolymer regions.
#' @param annotated_region_only logical(1): whether to only consider variants outside
#' annotated regions. If \code{TRUE}, only variants outside annotated regions will be
#' returned. If \code{FALSE}, all variants will be returned, which could take significantly
#' longer time.
#' @param names_from character(1): the column name in the metadata column of the annotation
#' (\code{mcols(annotation)[, names_from]}) to use for the \code{region} column in the output.
#' @return A tibble with columns: seqnames, pos, nucleotide, count, sum, freq, ref, region,
#' homopolymer_pct, bam_path The homopolymer percentage is calculated as the percentage of the
#' most frequent nucleotide in a window of \code{homopolymer_window} nucleotides around
#' the variant position, excluding the variant position itself.
#' @examples
#' ppl <- example_pipeline("SingleCellPipeline")
#' ppl <- run_step(ppl, "genome_alignment")
#' variants <- find_variants(
#'   bam_path = ppl@genome_bam,
#'   reference = ppl@genome_fa,
#'   annotation = ppl@annotation,
#'   min_nucleotide_depth = 4
#' )
#' head(variants)
#' @export
find_variants <- function(bam_path, reference, annotation, min_nucleotide_depth = 100,
                          homopolymer_window = 3, annotated_region_only = FALSE,
                          names_from = "gene_name", threads = 1) {
  if (is.character(reference)) {
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Reading reference ..."))
    reference <- Biostrings::readDNAStringSet(reference)
    # get rid of the `1` from >chr1 1
    names(reference) <- sapply(names(reference), function(x) strsplit(x, " ")[[1]][1])
  }
  if (is.character(annotation)) {
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Reading annotation ..."))
    annotation <- rtracklayer::import(annotation) |>
      (\(x) x[S4Vectors::mcols(x)$type == "gene", ])()
  }

  message("Merging overlapping genes ...")
  S4Vectors::mcols(annotation) <- S4Vectors::mcols(annotation)[, names_from, drop = FALSE]
  disjoint_annotation <- GenomicRanges::disjoin(annotation,
    ignore.strand=TRUE, with.revmap=TRUE
  )
  merge_mcol <- function(revmap, mcol) {
    sapply(revmap, function(x) {
      if (length(x) == 1) {
        return(mcol[x])
      } else {
        return(paste0(mcol[x], collapse = ", "))
      }
    })
  }
  message(sprintf("%d overlapping regions, their %s(s) are merged with `, ` as separator",
    sum(sapply(S4Vectors::mcols(disjoint_annotation)$revmap, length) > 1), names_from
  ))
  S4Vectors::mcols(disjoint_annotation) <- S4Vectors::mcols(disjoint_annotation)$revmap |>
    merge_mcol(S4Vectors::mcols(annotation)[, names_from]) |>
    S4Vectors::DataFrame() |>
    setNames(names_from)
  annotation <- disjoint_annotation

  if (!annotated_region_only) {
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Adding unannotated gaps ..."))
    # parsed annotation might not have seqlengths in seqinfo
    if (!all(as.character(GenomeInfoDb::seqnames(annotation)) %in% names(reference))) {
      warning("Some seqnames in annotation not found in reference")
      annotation <- annotation[as.character(GenomeInfoDb::seqnames(annotation)) %in% names(reference)]
    }
    GenomeInfoDb::seqinfo(annotation) <-
      reference[GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(annotation))] |>
      Biostrings::seqinfo()
    if (any(is.na(GenomeInfoDb::seqlengths(annotation)))) {
      stop("Missing seqlengths in seqinfo of annotation, cannot add unannotated gaps")
    }
    annotation <- c(annotation, GenomicRanges::gaps(annotation, ignore.strand = TRUE))
  }

  if (length(bam_path) == 1) {
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Got 1 bam file, parallelizing over each region ..."))
    variants <- tryCatch({
      BiocParallel::bplapply(
        sapply(seq_along(annotation), function(x) annotation[x]),
        function(grange) {
          find_variants_grange(bam_path, reference, grange, min_nucleotide_depth, names_from)
        },
        BPPARAM = BiocParallel::MulticoreParam(
          workers = threads, stop.on.error = TRUE, progressbar = TRUE # interactive()
        )
      )}, error = identity)
  } else {
    # multi-threading over bam files x granges
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Got multiple bam files, preparing for multi-threading ..."))
    args_grid <- expand.grid(
      grange = seq_along(annotation), # does not work on GRanges directly
      bam = bam_path,
      stringsAsFactors = FALSE
    )

    message(paste0(format(Sys.time(), "%H:%M:%S "), "Multi-threading over bam files x ranges ..."))
    variants <- tryCatch({
      BiocParallel::bpmapply(
        function(bam, grange) {
          find_variants_grange(bam, reference, annotation[grange], min_nucleotide_depth, names_from)
        },
        bam = args_grid$bam, grange = args_grid$grange,
        SIMPLIFY = FALSE,
        BPPARAM = BiocParallel::MulticoreParam(
          workers = threads, stop.on.error = TRUE, progressbar = TRUE # interactive()
        )
      )
    }, error = identity)
  }

  if (inherits(variants, "error")) {
    warning("Error occurred in `find_variants`, returning error object")
    return(variants)
  }

  message(paste0(format(Sys.time(), "%H:%M:%S "), "Merging results ..."))
  variants <- dplyr::bind_rows(variants)

  if (homopolymer_window > 1) {
    if (nrow(variants) == 0) {
      variants <- variants |>
        dplyr::mutate(homopolymer_pct = numeric(0))
    } else {
      message(paste0(format(Sys.time(), "%H:%M:%S "), "Calculating homopolymer percentages ..."))
      variants$homopolymer_pct <- homopolymer_pct(
        reference, variants$seqnames, variants$pos,
        include_alt = FALSE, n = homopolymer_window, threads = threads
      )
    }
  }

  return(variants)
}

#' mutation positions within the gene body
#'
#' Given a set of mutations and a gene annotation, calculate the position of each mutation
#' within the gene body. The gene annotation must have the following types: "gene" and "exon".
#' The gene annotation must be for one gene only. The mutations must be within the gene region.
#' The function will merge overlapping exons and calculate the position of each mutation
#' within the gene body, excluding intronic regions.
#'
#' @keywords internal
#' @importFrom GenomicRanges GRanges findOverlaps reduce
#' @importFrom IRanges IRanges start end width
#' @importFrom S4Vectors mcols queryHits subjectHits
#' @importFrom magrittr equals
#' @importFrom BiocGenerics strand
#'
#' @param mutations either the tibble output from \code{find_variants} or a GRanges object.
#' Make sure to filter it for only the gene of interest.
#' @param annotation_grange GRanges: the gene annotation. Must have the following types: "gene" and "exon".
#' @param verbose logical(1): whether to print messages.
#' @param type character(1): the type of position to calculate. Can be one of "TSS" (distance from the
#' transcription start site), "TES" (distance from the transcription end site), or "relative" (relative
#' position within the gene body).
#' @return A numeric vector of positions of each mutation within the gene body. When \code{type = "relative"},
#' the positions are normalized to the gene length, ranging from 0 (start of the gene) to 1 (end of the gene).
#' When \code{type = "TSS"} / \code{type = "TES"}, the distances from the transcription start
#' / end site.
mutation_positions_single <- function(mutations, annotation_grange, type, verbose = TRUE) {

  if (is.data.frame(mutations)) {
    # verify that all mutations are in the same gene
    if (!length(unique(mutations$region)) == 1) {
      errMsg <- sprintf(
        "Incorrent number of unique values in `mutations` for column `region`: %s",
        paste0(unique(mutations$region), collapse = ", ")
      )
      stop(errMsg)
  }
   
    # convert to GRanges of mutations 
    mutations <- GenomicRanges::GRanges(
      seqnames = mutations$seqnames,
      ranges = IRanges::IRanges(start = mutations$pos, end = mutations$pos)
    )
  } else {
    stopifnot("muations must be either a GRanges or a data.frame object" = is(mutations, "GRanges"))
  }

  # verify that the gene annotation is for one gene only
  if (!length(unique(S4Vectors::mcols(annotation_grange)$gene_id)) == 1) {
    errMsg <- sprintf(
      "Incorrent number of gene_id(s) in `annotation_grange`: %s",
      paste0(unique(S4Vectors::mcols(annotation_grange)$gene_id), collapse = ", ")
    )
    stop(errMsg)
  }

  # verify that all mutations are within the gene region
  mutations_ok <- GenomicRanges::findOverlaps(
    mutations,
    subset(annotation_grange, type == "gene")
  ) |>
    S4Vectors::queryHits() |>
    length() |>
    magrittr::equals(length(mutations))
  if (!mutations_ok) {
    stop("Some mutations are not within the gene region")
  }

  # merge overlapping exons
  merged_exons <- annotation_grange |>
    subset(type == "exon") |>
    GenomicRanges::reduce()

  overlaps <- GenomicRanges::findOverlaps(
    mutations, merged_exons,
    type = "within", select = "first"
  )
  mutations_in_exons <- mutations[!is.na(overlaps)]
  overlapped_idx <- overlaps[!is.na(overlaps)]
  if (verbose) {
    message(
      sprintf(
        "Found %d mutations in exons out of at total of %d mutations in %s",
        length(mutations_in_exons), length(mutations),
        unique(S4Vectors::mcols(annotation_grange)$gene_name)
      )
    )
  }
  # no mutations in exons
  if (length(mutations_in_exons) == 0) {
    return(rep(NA, length(mutations)))
  }

  cumulative_exon_lengths <- cumsum(IRanges::width(merged_exons))

  positions <-
    mapply(
      function(mutation_id, exon_idx) {
        mutation <- mutations_in_exons[mutation_id]
        exon_end <- IRanges::end(merged_exons[exon_idx])
        exon_cumulative_length <- cumulative_exon_lengths[exon_idx]
        position_tss <- exon_cumulative_length - exon_end + IRanges::start(mutation)
        if (type %in% c("TSS", "TES")) {
          return(position_tss)
        } else if (type == "relative") {
          return(position_tss / sum(IRanges::width(merged_exons)))
        }
      },
      mutation_id = seq_along(mutations_in_exons), exon_idx = overlapped_idx,
      SIMPLIFY = TRUE
    )

  # flip the positions if the gene is on the negative strand
  neg_strand <- subset(annotation_grange, type == "gene") |>
    BiocGenerics::strand() |>
    as.character() |>
    magrittr::equals("-")
  if (neg_strand && type == "relative") {
    positions <- 1 - positions
  } else if ((neg_strand && type == "TSS") || (!neg_strand && type == "TES")) {
    positions <- sum(IRanges::width(merged_exons)) - positions
  }

  # map the positions back to the original mutations (and add NAs)
  full_positions <- rep(NA, length(mutations))
  full_positions[!is.na(overlaps)] <- positions

  return(full_positions)
}

#' Calculate mutation positions within the gene body
#'
#' Given a set of mutations and gene annotation, calculate the position of each mutation
#' within the gene body they are in.
#'
#' @param mutations either the tibble output from \code{find_variants}. It must have columns \code{seqnames},
#' \code{pos}, and a third column for specifying the gene id or gene name. The mutation must be within the gene region.
#' @param annotation Either path to the annotation file (GTF/GFF) or a GRanges object of the gene annotation.
#' @param type character(1): the type of position to calculate. Can be one of "TSS" (distance from the
#' transcription start site), "TES" (distance from the transcription end site), or "relative" (relative
#' position within the gene body).
#' @param bin logical(1): whether to bin the relative positions into 100 bins. Only applicable when
#' \code{type = "relative"}.
#' @param by character(1): the column name in the annotation to match with the gene annotation.
#' E.g. \code{c("region" = "gene_name")} to match the `region` column in the mutations with the
#' `gene_name` column in the annotation.
#' @param threads integer(1): number of threads to use.
#' @return A numeric vector of positions of each mutation within the gene body. When \code{type = "relative"},
#' the positions are normalized to the gene length, ranging from 0 (start of the gene) to 1 (end of the gene).
#' When \code{type = "TSS"} / \code{type = "TES"}, the distances from the transcription start
#' / end site. If \code{bin = TRUE}, and \code{type = "relative"}, the relative positions are binned into 100
#' bins along the gene body, and the output is a matrix with the number of mutations in each bin, the
#' rows are named by the \code{by} column (e.g. gene name).
#' @examples
#' variants <- data.frame(
#'   seqnames = rep("chr14", 8),
#'   pos = c(1084, 1085, 1217, 1384, 2724, 2789, 5083, 5147),
#'   region = rep("Rps24", 8)
#' )
#' positions <-
#'  mutation_positions(
#'    mutations = variants,
#'    annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
#'  )
#' @importFrom rtracklayer import
#' @importFrom dplyr mutate filter
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom S4Vectors mcols
#' @export
mutation_positions <- function(mutations, annotation, type = "relative", bin = FALSE, by = c("region" = "gene_name"), threads = 1) {
  if (!length(by) == 1) {
    stop("'by' must be a character vector of length 1")
  }
  if (is.null(names(by))) {
    names(by) <- by
  }

  stopifnot(
    "Some mutations do not have a value for the 'by' column" = !any(is.na(mutations[, names(by)]))
  )

  stopifnot(
    "'type' must be one of 'TSS', 'TES', 'relative'" = type %in% c("TSS", "TES", "relative")
  )

  if (is.character(annotation)) {
    message(paste0(format(Sys.time(), "%H:%M:%S "), "Reading annotation ..."))
    annotation_grange <- rtracklayer::import(annotation)
  } else {
    annotation_grange <- annotation
  }

  stopifnot(
    "Some gene(s) are not found in the annotation" =
      all(mutations[, names(by), drop = TRUE] %in% S4Vectors::mcols(annotation_grange)[, by])
  )

  mutations_split <- mutations |>
    dplyr::mutate(row_id = row_number(),
      region = mutations[, names(by), drop = TRUE]) |>
    split(mutations[, names(by), drop = TRUE])
  # idx to map output back to original order of mutations
  mutations_split_idx <- lapply(mutations_split, function(x) x$row_id) |>
    unlist() |>
    match(x = seq_len(nrow(mutations)), table = _)

  coverages <- tryCatch({
    BiocParallel::bplapply(
      names(mutations_split),
      function(x) {
        annot_i <- subset(annotation_grange, S4Vectors::mcols(annotation_grange)[, by] == x)
        pos <- mutation_positions_single(mutations_split[[x]], annot_i, type = type, verbose = FALSE)
        if (!bin || type != "relative") {
          return(pos)
        }
        pos <- round(pos * 100)
        coverage <- sapply(1:100, function(i) sum(pos == i, na.rm = TRUE))
        return(coverage)
      },
      BPPARAM = BiocParallel::MulticoreParam(
        workers = threads, stop.on.error = TRUE, progressbar = TRUE
      )
    )
  }, error = identity)
  if (inherits(coverages, "error")) {
    warning("Error occurred in `mutations_coverage`, returning error object")
    return(coverages)
  }

  if (!bin) {
    names(coverages) <- names(mutations_split)
    coverages <- unlist(coverages)
    coverages <- coverages[mutations_split_idx]
    return(coverages)
  }

  # sum up all coverages
  total_coverage <- do.call(rbind, coverages)
  rownames(total_coverage) <- names(mutations_split)
  return(total_coverage)
}
