#!/bin/bash

#Download cutadapt  

## make sure bioconda is set up 
### setup requirements have changed in the last few years, so this should be done regardless of when it was downloaded
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

## update or install cutadapt from bioconda
conda create -n cutadapt cutadapt
## check version 
conda activate cutadapt
cutadapt --version 
# 5.0

# set directories 
## $HOME = home dir
FILES="$HOME/fastq"
OUT="$HOME/trimmed"
cutadapt="/Users/lauragivens/miniconda3/envs/cutadapt"
# primer sequences
#citation:
FwdPrimer=""
RevPrimer=""

