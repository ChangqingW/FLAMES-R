FROM bioconductor/bioconductor_docker:devel

WORKDIR /home/rstudio

# Install system dependencies first
RUN apt-get update && apt-get install -y \
    samtools \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy source files
COPY --chown=rstudio:rstudio . /home/rstudio/

RUN Rscript -e "options(Ncpus = 8); BiocManager::install(ask=FALSE); devtools::install('.', dependencies=TRUE, build_vignettes=TRUE, repos = BiocManager::repositories())"

# Pre-initialize basilisk environment
RUN su - rstudio -c "Rscript -e 'basilisk::basiliskRun(env = FLAMES:::flames_env, fun = function(){})'" \
    && chmod -R 777 /home/rstudio/.cache/R/basilisk
