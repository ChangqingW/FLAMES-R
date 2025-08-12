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
RUN Rscript -e "install.packages('remotes'); remotes::install_deps('.', dependencies=TRUE)"

# copy everything else
COPY . .
# install FLAMES package
RUN Rscript -e "BiocManager::install('basilisk', type = 'source', force = TRUE); BiocManager::install('GoekeLab/bambu', force = TRUE); remotes::install_local('.', dependencies=TRUE)"

# test to see if FLAMES is installed
# if the import succeeds, this command will exit with code 0
RUN Rscript -e "library(FLAMES)" && \
    install_path=$(Rscript -e 'cat(system.file(package = "FLAMES"))') && \
    mkdir -p "$install_path/bin" && \
    wget -qO - https://github.com/COMBINE-lab/oarfish/releases/download/v0.8.1/oarfish-x86_64-unknown-linux-gnu.tar.xz | \
    tar -xJ --strip-components=1 -C "$install_path/bin" 'oarfish-x86_64-unknown-linux-gnu/oarfish'

RUN Rscript -e 'basilisk::basiliskRun(env = FLAMES:::flames_env, fun = function(){})'

CMD ["R", "--no-save"]
