#!/bin/bash
conda activate fastqc 

#move report files from cutadapt into their own directory
cd /Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI
mkdir -p cutadapt/report
 mv fastq/*.json cutadapt/report/
#move undetermined barcode files into their own directory
mkdir -p UnknownBarcodes/trim
 mv cutadapt/trim/Undetermined* UnknownBarcodes/trim
 
 #go to 02_Trim_QC.R
 #plotqualityprofile for R reads fails with BiocParallel error; Error in density.default(qscore):'x' contains missing values
 #my guess is this is because of the reports from cutadapt showing that such a low percentage of r reads are passing 
 
#fastqc and multiqc report
 RAWDIR=/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI
 TRIMDIR=$RAWDIR/cutadapt/FQCreports
 mkdir -p $TRIMDIR
 
 find $RAWDIR/cutadapt/trim -type f -name "*.fastq" -o -name "*.fastq.gz" | xargs fastqc --threads 4 outdir "$TRIMDIR"

#aggregate all the fastqc files into a summary file 
cd $TRIMDIR

multiqc . #run multiqc

#1 Apr 2025

