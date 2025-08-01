% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-BulkPipeline.R
\name{steps<-}
\alias{steps<-}
\alias{steps<-,FLAMES.Pipeline-method}
\title{Set steps to perform in the pipeline}
\usage{
steps(pipeline) <- value

\S4method{steps}{FLAMES.Pipeline}(pipeline) <- value
}
\arguments{
\item{pipeline}{An object of class `FLAMES.Pipeline`}

\item{value}{A named logical vector containing all possible steps
for the pipeline. The names of the vector are the step names,
and the values are logical indicating whether the step is
configured to be performed.}
}
\value{
An pipeline of class `FLAMES.Pipeline` with the updated steps.
}
\description{
Set steps to perform in the pipeline
}
\examples{
ppl <- example_pipeline()
steps(ppl) <- c(
  barcode_demultiplex = TRUE,
  genome_alignment = TRUE,
  gene_quantification = TRUE,
  isoform_identification = FALSE,
  read_realignment = FALSE,
  transcript_quantification = TRUE
)
ppl
# or partially change a step:
steps(ppl)["read_realignment"] <- TRUE
ppl
}
