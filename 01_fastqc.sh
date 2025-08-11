#!/bin/bash

############################# Install #############################  

# FastQC  https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt
conda create --name fastqc bioconda::fastqc
# activate the environment  
conda activate fastqc 
# install multiqc inside the environment
# conda install multiqc #does not seem to work with m1 mac
conda install -c bioconda -c conda-forge multiqc #for M1 Mac
pip install importlib-metadata
# should now have both FastQC and MultiQC in the same conda environment  

#check versions  
fastqc --version 
#v0.12.1
multiqc --version
#v1.10.1


# First, we are going to look at the quality of our reads  
RAWDIR=/Sequences #set the RAWDIR variable to the location of the sequencing files  
OUTDIR=$RAWDIR/FQCreports_raw #name the folder where you want the fastqc reports to go  
mkdir -p $OUTDIR #make the OUTDIR folder

find $RAWDIR/fastq -type f -name "*.fastq" -o -name "*.fastq.gz" | xargs fastqc --threads 4 --outdir "$OUTDIR" 
#makes two files for each fastq.gz file; an html and a zip file

#aggregate all the fastqc files into a summary file 
cd $OUTDIR
mv Undetermined_* ../UnknownBarcodes #remove undetermined reads from multiqc analysis

multiqc . #run multiqc
mv multiqc* ../ #move the multiqc files out of the FQCreports file 

#1 Apr 2025

