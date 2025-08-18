# Multi-stage build for better caching
# Stage 1: Build base with all dependencies
FROM rocker/r-ver:4.5.1 AS deps-builder

WORKDIR /tmp/build

# Install system dependencies for building R packages
RUN apt-get update && \
    apt-get install -y \
        samtools \
        libglpk-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libpng-dev \
        libmagick++-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        zlib1g-dev \
        pkg-config \
        wget \
        git \
        ca-certificates \
        build-essential \
        gfortran && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set Bioconductor version for devel
ENV BIOCONDUCTOR_VERSION=3.22
ENV BASILISK_USE_SYSTEM_DIR=1

# Install BiocManager and set to devel version
RUN Rscript -e 'install.packages("BiocManager", repos="https://cran.rstudio.com")' && \
    Rscript -e "BiocManager::install(version='${BIOCONDUCTOR_VERSION}', update=TRUE, ask=FALSE)" && \
    Rscript -e "install.packages('remotes', repos=BiocManager::repositories())"

# Copy only DESCRIPTION for dependency installation
COPY DESCRIPTION .

# Install all dependencies (this creates a reusable layer)
RUN Rscript -e "remotes::install_deps('.', dependencies=TRUE, repos=BiocManager::repositories())"

# Stage 2: Final image
FROM rocker/r-ver:4.5.1

WORKDIR /home/flames

# Install runtime dependencies (based on working container analysis)
RUN apt-get update && \
    apt-get install -y \
        samtools \
        libglpk-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libpng-dev \
        libmagick++-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        zlib1g-dev \
        wget \
        git \
        ca-certificates && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy R libraries from builder stage
COPY --from=deps-builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library

# Set environment variables
ENV BIOCONDUCTOR_VERSION=3.22
ENV BASILISK_USE_SYSTEM_DIR=1

# Copy project files
COPY . .

# Install only FLAMES package (dependencies already installed)
RUN Rscript -e "BiocManager::install('basilisk', type = 'source', force = TRUE); install.packages('.', repos = NULL, type = 'source')"

# Test installation and download oarfish
RUN Rscript -e "library(FLAMES)" && \
    install_path=$(Rscript -e 'cat(system.file(package = "FLAMES"))') && \
    mkdir -p "$install_path/bin" && \
    wget -qO - https://github.com/COMBINE-lab/oarfish/releases/download/v0.8.1/oarfish-x86_64-unknown-linux-gnu.tar.xz | \
    tar -xJ --strip-components=1 -C "$install_path/bin" 'oarfish-x86_64-unknown-linux-gnu/oarfish'

# Initialize basilisk environment
RUN Rscript -e 'basilisk::basiliskRun(env = FLAMES:::flames_env, fun = function(){})'

CMD ["R", "--no-save"]
