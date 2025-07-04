#' Isoform identification
#' @description Long-read isoform identification with FLAMES or bambu.
#' @param annotation Path to annotation file. If configured to use bambu, the annotation
#' must be provided as GTF file.
#' @param genome_fa The file path to genome fasta file.
#' @param genome_bam File path to BAM alignment file. Multiple files could be provided.
#' @param outdir The path to directory to store all output files.
#' @param config Parsed FLAMES configurations.
#' @return The updated annotation and the transcriptome assembly will be saved in the
#' output folder as \code{isoform_annotated.gff3} (GTF if bambu is selected) and
#' \code{transcript_assembly.fa} respectively.
#' @keywords internal
find_isoform <- function(annotation, genome_fa, genome_bam, outdir, config) {
  # pipeline types: singe_cell, single_cell_multisample, bulk
  if (config$pipeline_parameters$bambu_isoform_identification) {
    find_isoform_bambu(annotation, genome_fa, genome_bam, outdir, config)
  } else {
    find_isoform_flames(annotation, genome_fa, genome_bam, outdir, config)
  }
}

#' @importFrom bambu writeToGTF prepareAnnotations bambu
#' @importFrom withr with_package
#' @importFrom SummarizedExperiment assays rowRanges
find_isoform_bambu <- function(annotation, genome_fa, genome_bam, outdir, config) {
  # if annotation is .gtf.gz, unzip as a temp file
  useTempAnnot <- FALSE
  bambuTempAnnot <- ""
  if (stringr::str_ends(annotation, ".gz")) {
    cat("Unzipping annotation file for bambu\n")
    useTempAnnot <- TRUE
    bambuTempAnnot <- R.utils::gunzip(annotation, remove = FALSE)
    annotation <- bambuTempAnnot # override using zipped annotation file
  }
  bambuAnnotations <- bambu::prepareAnnotations(annotation)
  # Tmp fix: remove withr if bambu imports seqlengths properly
  # https://github.com/GoekeLab/bambu/issues/255
  # min.readCount seems to cause errors
  # https://github.com/GoekeLab/bambu/issues/364

  bambu_out <- withr::with_package("GenomeInfoDb",
    bambu::bambu(
      reads = genome_bam,
      annotations = bambuAnnotations,
      genome = genome_fa,
      quant = TRUE,
      discovery = ifelse(
        is.null(config$isoform_parameters$bambu_discovery),
        TRUE,
        config$isoform_parameters$bambu_discovery),
      lowMemory = TRUE,
      NDR = config$isoform_parameters$bambu_ndr,
      verbose = config$isoform_parameters$bambu_verbose,
      # https://github.com/GoekeLab/bambu/issues/416#issuecomment-1987499886
      # ncore = config$pipeline_parameters$threads
  ))

  bambu::writeToGTF(SummarizedExperiment::rowRanges(bambu_out), file.path(outdir, "isoform_annotated_unfiltered.gtf"))
  if (is.null(config$isoform_parameters$bambu_trust_reference) || config$isoform_parameters$bambu_trust_reference) {
    bambu_out <- bambu_out[base::rowSums(SummarizedExperiment::assays(bambu_out)$counts) >= 1, ]
  } else {
    bambu_out <- bambu_out[base::rowSums(SummarizedExperiment::assays(bambu_out)$counts) >= config$isoform_parameters$min_sup_cn, ]
  }

  isoform_gtf <- file.path(outdir, "isoform_annotated.gtf") # Bambu outputs GTF
  bambu::writeToGTF(SummarizedExperiment::rowRanges(bambu_out), isoform_gtf)
  # annotation_to_fasta(isoform_gtf, genome_fa, outdir)

  if (useTempAnnot & file.exists(bambuTempAnnot)) {
    file.remove(bambuTempAnnot)
  }

  return(isoform_gtf)
}

#' @importFrom reticulate import_from_path
#' @importFrom Rsamtools indexFa
#' @importFrom basilisk basiliskRun
find_isoform_flames <- function(annotation, genome_fa, genome_bam, outdir, config) {
  if (!all(file.exists(paste0(genome_bam, ".bai")))) {
    stop(c("Cannot find corresponding bam file(s) ", paste0(genome_bam, ".bai"), ". Cancelling find_isoform."))
  }
  isoform_annotation <- file.path(outdir, "isoform_annotated.gff3")
  tss_stat <- file.path(outdir, "tss_tes.bedgraph")
  transcript_assembly <- file.path(outdir, "transcript_assembly.fa")
  raw_splice <- file.path(outdir, "splice_raw.gff3")

  if (length(genome_bam) == 1) {
    if (config$pipeline_parameters$multithread_isoform_identification) {
      # C++ Multithreaded implementation of python find_isoform
      find_isoform_multithread(
        annotation, genome_bam, isoform_annotation, tss_stat, genome_fa, transcript_assembly, config$isoform_parameters, ifelse(config$isoform_parameters$generate_raw_isoform, raw_splice, "")
      )
      annotation_to_fasta(isoform_annotation, genome_fa, outdir)
    } else {
      ret <- basiliskRun(
        env = flames_env, fun = function(gff3, genome, iso, tss, fa, tran, ds, conf, raw) {
            python_path <- system.file("python", package = "FLAMES")
            find <- reticulate::import_from_path("find_isoform", python_path)
            ret <- find$find_isoform(gff3, genome, iso, tss, fa, tran, ds, conf, raw)
            ret
        },
        gff3 = annotation,
        genome = genome_bam,
        iso = isoform_annotation,
        tss = file.path(outdir, "tss_tes.bedgraph"),
        fa = genome_fa,
        tran = transcript_assembly,
        ds = config$isoform_parameters$downsample_ratio,
        conf = config, 
        raw = ifelse(config$isoform_parameters$generate_raw_isoform, file.path(outdir, "splice_raw.gff3"), FALSE)
      )
    }
  } else {
    ret <- basiliskRun(
        env = flames_env,
        fun = function(gff3, genome, iso, tss, fa, tran, ds, conf, raw) {
          python_path <- system.file("python", package = "FLAMES")
          find <- reticulate::import_from_path("find_isoform", python_path)
          ret <- find$find_isoform_multisample(gff3, genome, iso, tss, fa, tran, ds, conf, raw)
          ret
        },
        gff3 = annotation,
        genome = genome_bam,
        iso = isoform_annotation,
        tss = tss_stat,
        fa = genome_fa,
        tran = transcript_assembly,
        ds = config$isoform_parameters$downsample_ratio,
        conf = config,
        raw = ifelse(config$isoform_parameters$generate_raw_isoform, raw_splice, FALSE)
    )
  }

  return(isoform_annotation)
}


#' Fake stranded GFF file
#' @description Check if all the transcript in the annotation is stranded. If not, convert to '+'.
#' @return Path to the temporary file with unstranded transcripts converted to '+'.
#' @keywords internal
fake_stranded_gff <- function(gff_file) {
  # check if all the transcript in the annotation is stranded
  annotation_d <- read.csv(gff_file, sep = "\t",
    header = FALSE, stringsAsFactors = FALSE,
    comment.char = "#")
  strands <- annotation_d[, 7]
  if (any(strands == '.')) {
    modified_gtf <- paste0(tempfile(), '/tmp.gtf')
    dir.create(dirname(modified_gtf))
    warning(sprintf("Some transcripts in the annotation file %s are not stranded. Converting to '+' in temporary file %s", gff_file, modified_gtf))
    strands[strands == '.'] <- '+'
    annotation_d[, 7] <- strands
    write.table(annotation_d, modified_gtf, sep = "\t",
      row.names = FALSE, quote = FALSE, col.names = FALSE)
    return(modified_gtf)
    # file will get deleted after quitting R
  }
  return(gff_file)
}

#' GTF/GFF to FASTA conversion
#' @description convert the transcript annotation to transcriptome assembly as FASTA file. The
#' genome annotation is first imported as TxDb object and then used to extract transcript sequence
#' from the genome assembly.
#' @param isoform_annotation Path to the annotation file (GTF/GFF3)
#' @param genome_fa The file path to genome fasta file.
#' @param outfile The file path to the output FASTA file.
#' @param extract_fn (optional) Function to extract \code{GRangesList} from the genome TxDb object.
#' E.g. \code{function(txdb){GenomicFeatures::cdsBy(txdb, by="tx", use.names=TRUE)}}
#' @return This does not return anything. A FASTA file will be created at the specified location.
#'
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom GenomicFeatures extractTranscriptSeqs
#' @importFrom txdbmaker makeTxDbFromGFF
#' @importFrom Rsamtools indexFa
#' @importFrom utils write.table
#'
#' @examples
#' fasta <- tempfile()
#' annotation_to_fasta(system.file("extdata", "rps24.gtf.gz", package = "FLAMES"), system.file("extdata", "rps24.fa.gz", package = "FLAMES"), fasta)
#' cat(readChar(fasta, 1e3))
#'
#' @export
annotation_to_fasta <- function(isoform_annotation, genome_fa, outfile, extract_fn) {
  isofile <- fake_stranded_gff(isoform_annotation)

  dna_string_set <- Biostrings::readDNAStringSet(genome_fa)
  names(dna_string_set) <- gsub(" .*$", "", names(dna_string_set))
  
  txdb <- txdbmaker::makeTxDbFromGFF(isofile)
  
  # Close the database connection on exit
  on.exit(DBI::dbDisconnect(txdb$conn), add = TRUE)

  if (missing(extract_fn)) {
    tr_string_set <- GenomicFeatures::extractTranscriptSeqs(dna_string_set, txdb,
      use.names = TRUE)
  } else {
    extracted_grl <- extract_fn(txdb)
    tr_string_set <- GenomicFeatures::extractTranscriptSeqs(dna_string_set, extracted_grl)
    # additional arguments are allowed only when 'transcripts' is not a GRangesList object
  }

  if (length(names(tr_string_set)) > length(unique(names(tr_string_set)))) {
    cat("Duplicated transcript IDs present, removing ...")
    tr_string_set <- tr_string_set[unique(names(tr_string_set))]
  }

  Biostrings::writeXStringSet(tr_string_set, outfile)
  Rsamtools::indexFa(outfile)
  return(invisible())
}

#' Parse FLAMES' GFF output
#' @description Parse FLAMES' GFF ouputs into a Genomic Ranges List
#' @param file the GFF file to parse
#' @return A Genomic Ranges List
#' @keywords internal
get_GRangesList <- function(file) {
  if (is.character(file)) {
    isoform_gr <- rtracklayer::import(file, feature.type = c("exon", "utr"))
    if (grepl("\\.gff3(\\.gz)?$", file)) {
      isoform_gr$Parent <- as.character(isoform_gr$Parent)
      isoform_gr$transcript_id <- unlist(lapply(strsplit(isoform_gr$Parent, split = ":"), function(x) {
        x[2]
      }))
    }
  } else if (is(file, "GRanges")) {
    isoform_gr <- file
  } else if (is(file, "GRangesList")) {
    return(file)
  } else {
    stop(sprintf("Unsupported input type: %s", class(file)))
  }

  isoform_grl <- S4Vectors::split(isoform_gr, isoform_gr$transcript_id)
  return(isoform_grl)
}
