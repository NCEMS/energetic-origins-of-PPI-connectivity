#!/bin/bash

conda create --name snakemake -c bioconda -y snakemake
#conda create --name centrality-cos-dist --file=centrality-cos-dist.yml
conda env create --name network-analysis --file=network-analysis.yml
conda env create --name cagiada-stability --file=cagiada-stability.yml
#conda env create --name deeptmhmm --file=deeptmhmm.yml
