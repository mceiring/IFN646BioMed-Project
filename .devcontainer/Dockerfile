# Load the images
FROM nanozoo/bowtie2
FROM r-base

# Install the required R packages (and dependencies)
RUN apt-get update && apt-get install -y libcurl4-openssl-dev
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev

RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager');\
BiocManager::install('polyester');\
BiocManager::install('Biostrings');\
BiocManager::install('rtracklayer');\
install.packages('jsonlite', repos='http://cran.rstudio.com/');\
install.packages('languageserver');"