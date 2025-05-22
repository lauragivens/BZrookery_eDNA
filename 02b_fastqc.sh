#!/bin/bash
conda activate fastqc 

#move report files from cutadapt into their own directory
cd /Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI
 
#fastqc and multiqc report
 DIR=/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI/cutadapt/trim
 REPORT=$DIR/reports-fastqc
 mkdir -p $REPORT
 
 find $DIR -type f -name "*.fastq" -o -name "*.fastq.gz" | xargs fastqc --threads 4 outdir "$REPORT"

#aggregate all the fastqc files into a summary file 
multiqc $DIR/reports-cutadapt/ --filename "multiqc-cutadapt-report" #

multiqc $REPORT/ --filename "multiqc-cutadapt-reads-report" #check overall quality of reads
multiqc $REPORT/*R1* --filename "$REPORT/multiqc-cutadapt-reads-report-R1" #check quality of forward reads
multiqc $REPORT/*R2* --filename "$REPORT/multiqc-cutadapt-reads-report-R2" #check quality of reverse reads

#both F and R start dropping off at 270 bp

#1 Apr 2025
#edited 2 Apr 2025
