setClass(
  "FLAMES.Pipeline",
  slots = list(
    # inputs
    config = "list",            # Configuration parameters
    outdir = "character",       # Output directory
    fastq = "character",        # Path to the FASTQ files
    annotation = "character",   # Path to the annotation file
    genome_fa = "character",    # Path to the genome FASTA file

    # outputs
    metadata = "list",          # Metadata for the pipeline
    genome_bam = "character",     # Path to the genome BAM file
    transcriptome_bam = "character", # Path to the transcript BAM file
    novel_isoform_annotation = "character", # Path to the novel isoform GFF / GTF file
    transcriptome_assembly = "character", # Path to the transcriptome assembly file
    experiment = "SummarizedExperiment", # SummarizedExperiment object for quantification results

    # binaries
    minimap2 = "character",     # Path to the minimap2 binary
    k8 = "character",           # Path to the k8 binary
    samtools = "character",     # Path to the samtools binary

    # pipeline state
    steps = "logical",           # Steps to perform
    completed_steps = "logical", # Completed steps
    durations = "difftime",      # Durations of each step
    last_error = "list"
  ),
  prototype = list(
    config = list(),
    outdir = NA_character_,
    fastq = NA_character_,
    annotation = NA_character_,
    genome_fa = NA_character_,

    metadata = list(),
    genome_bam = NA_character_,
    transcriptome_bam = NA_character_,
    novel_isoform_annotation = NA_character_,
    transcriptome_assembly = NA_character_,
    # experiment = NA,

    minimap2 = NA_character_,
    k8 = NA_character_,
    samtools = NA_character_,

    steps = logical(),
    completed_steps = logical(),
    durations = structure(numeric(0), class = "difftime", units = "secs"),
    last_error = list()
  )
)

setClass(
  "FLAMES.SingleCellPipeline",
  contains = "FLAMES.Pipeline",
  slots = list(
    # inputs
    barcodes_file = "character",        # Path to the barcodes file
    expect_cell_number = "numeric",     # Expected number of cells
    # outputs
    demultiplexed_fastq = "character",  # path to demultiplexed FASTQ files
    deduped_fastq = "character"       # path to deduplicated FASTQ files
  ),
  prototype = list(
    barcodes_file = NA_character_,
    expect_cell_number = NA_real_,
    demultiplexed_fastq = NA_character_,
    deduped_fastq = NA_character_
  )
)

# Multi-sample pipeline is the same as SingleCellPipeline
# but input slots will be vectors of the same length
setClass(
  "FLAMES.MultiSampleSCPipeline",
  contains = "FLAMES.SingleCellPipeline",
  slots = list(
    experiments = "list" # List of SummarizedExperiment objects for each sample
  ),
  prototype = list(
    experiments = list()
  )
)
