% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-BulkPipeline.R
\name{steps}
\alias{steps}
\alias{steps,FLAMES.Pipeline-method}
\title{Steps to perform in the pipeline}
\usage{
steps(pipeline)

\S4method{steps}{FLAMES.Pipeline}(pipeline)
}
\arguments{
\item{pipeline}{An object of class `FLAMES.Pipeline`}
}
\value{
A named logical vector containing all possible steps
for the pipeline. The names of the vector are the step names,
and the values are logical indicating whether the step is
configured to be performed.
}
\description{
Steps to perform in the pipeline
}
\examples{
ppl <- example_pipeline()
steps(ppl)
}
