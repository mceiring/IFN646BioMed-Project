# Use Miniconda3 as the base image
FROM continuumio/miniconda3

# Create a new Conda environment with Python 3.8 and R
RUN conda create -n biomed python=3.8 r-base

# Activate the new environment and set it as the default for the shell
RUN echo "conda activate biomed" >> ~/.bashrc

# Refresh the shell
SHELL ["/bin/bash", "--login", "-c"]

# Install required packages, including bowtie2 and other bioinformatics tools
RUN conda install -c bioconda bowtie2
RUN conda install -c bioconda subread
RUN conda install -c bioconda kallisto

# Configure CRAN mirror for R
RUN echo 'options(repos = c(CRAN = "http://cran.rstudio.com/"))' >> ~/.Rprofile

# Install required R packages
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager');\
BiocManager::install('polyester');\
BiocManager::install('Biostrings');\
BiocManager::install('rtracklayer');\
BiocManager::install('ShortRead');\
BiocManager::install('DESeq2');\
install.packages('jsonlite');\
install.packages('fs', version='1.5.0');\
install.packages('pkgload', version='1.2.4');\
install.packages('roxygen2', version='7.1.2');\
install.packages('languageserver');"
# still have languageserver error (skip for now)