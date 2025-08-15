FROM rocker/r2u:24.04

WORKDIR /home/flames

# Install system dependencies first
RUN apt-get update && \
    apt-get install -y samtools libcurl4-openssl-dev && \
    apt-get autoremove -y && \
    apt-get clean

# install dependencies
ENV BASILISK_USE_SYSTEM_DIR=1
COPY DESCRIPTION .

# Need to remove the sed command after the next release when Seqinfo is available
RUN sed -i 's/Seqinfo/GenomeInfoDb/g' DESCRIPTION && \
    Rscript -e "install.packages('remotes'); remotes::install_deps('.', dependencies=TRUE)"
# RUN Rscript -e "install.packages('remotes'); remotes::install_deps('.', dependencies=TRUE)"

# copy everything else
COPY . .

# Remove this layer when Seqinfo is available in Bioconductor
RUN sed -i 's/Seqinfo/GenomeInfoDb/g' DESCRIPTION && \
    sed -i 's/Seqinfo/GenomeInfoDb/g' NAMESPACE && \
    find R/ -name "*.R" -type f -exec sed -i 's/Seqinfo/GenomeInfoDb/g' {} \;

# install FLAMES package
RUN Rscript -e "BiocManager::install('basilisk', type = 'source', force = TRUE); remotes::install_local('.', dependencies=TRUE)"

# test to see if FLAMES is installed
# if the import succeeds, this command will exit with code 0
RUN Rscript -e "library(FLAMES)" && \
    install_path=$(Rscript -e 'cat(system.file(package = "FLAMES"))') && \
    mkdir -p "$install_path/bin" && \
    wget -qO - https://github.com/COMBINE-lab/oarfish/releases/download/v0.8.1/oarfish-x86_64-unknown-linux-gnu.tar.xz | \
    tar -xJ --strip-components=1 -C "$install_path/bin" 'oarfish-x86_64-unknown-linux-gnu/oarfish'

RUN Rscript -e 'basilisk::basiliskRun(env = FLAMES:::flames_env, fun = function(){})'

CMD ["R", "--no-save"]
