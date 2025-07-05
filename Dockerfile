FROM rocker/r2u:24.04

WORKDIR /home/flames

# Install system dependencies first
RUN apt-get update && \
    apt-get install -y samtools rustup && \
    apt-get clean 
RUN rustup default stable
RUN Rscript -e "install.packages('remotes')"

# copy the R package description file so the dependency install can be cached
COPY DESCRIPTION .

# install dependencies
RUN Rscript -e "remotes::install_deps('.', dependencies=TRUE)" && apt-get autoremove -y && apt-get clean

# install dependencies for the FLAMES compilation
RUN apt-get install -y libcurl4-openssl-dev && \
    apt-get clean

# copy everything else
COPY . .
# install FLAMES package
RUN Rscript -e "remotes::install_local('.', dependencies=TRUE)"

# test to see if FLAMES is installed
# if the import succeeds, this command will exit with code 0
RUN Rscript -e "library(FLAMES)"

CMD ["R", "--no-save"]