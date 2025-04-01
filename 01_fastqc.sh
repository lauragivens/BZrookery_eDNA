#!/bin/bash
conda activate fastqc 
conda install multiqc
fastqc --version 
#v0.12.1

# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt
# 

# First, we are going to look at the quality of our reads  
RAWDIR=/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI
OUTDIR=$RAWDIR/FQCreports
mkdir -p $OUTDIR

find $RAWDIR/fastq -type f -name "*.fastq" -o -name "*.fastq.gz" | xargs fastqc --threads 4 --outdir "$OUTDIR" 
#makes two files for each fastq.gz file; an html and a zip file

#aggregate all the fastqc files into a summary file 
cd $OUTDIR
mv Undetermined_* ../UnknownBarcodes #remove undetermined reads from multiqc analysis

multiqc . #run multiqc
mv multiqc* ../ #move the multiqc files out of the FQCreports file 

#1 Apr 2025

