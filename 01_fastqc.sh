#!/bin/bash
conda activate fastqc 
fastqc --version 
#v0.12.1

# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt
# 
RAWDIR=/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI
OUTDIR=$RAWDIR/FQCreports
mkdir -p $OUTDIR

find $RAWDIR/fastq -type f -name "*.fastq" -o -name "*.fastq.gz" | xargs fastqc --threads 4 --outdir "$OUTDIR" 


#1 Apr 2025

