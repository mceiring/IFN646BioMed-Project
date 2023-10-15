#!/bin/bash

# Initialize conda
conda init bash

# Create a new conda environment and install packages
conda create -n bio_env python=3.9 -y
source activate bio_env
conda install -c bioconda bowtie2 salmon -y