#!/bin/bash

#Assign Taxonomy
BLASTDB="/Users/lauragivens/blastdb/core_nt"
export BLASTDB

dir_results=/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI/cutadapt/results
uniqueseqs="$dir_results"/dada2-uniqueseqs.fasta

blastn -db core_nt \
-query $uniqueseqs \
-perc_identity 97 \
-word_size 30 \
-outfmt 6 \
-num_threads 4 \
-out $dir_results/dada2.uniques.BLAST.default.txt

# this has been stopping on error 
# Error memory mapping:/Users/lauragivens/blastdb/core_nt/core_nt.17.nnd openedFilesCount=251 threadID=0
# BLAST Database error: Error pre-fetching sequence data 

# Also attempted with 
# changing number of threads from 4; tried 2 and 1 
# BLASTDB="/Volumes/Fuji/blastdb/core_nt" ## in case the issue was memory/RAM 
# and blastn -db MAR_taxid ## using the alias masked database of only MAR fish spp taxids 

## now transferring the core_nt database to DCC to see if it will work on a higher performing cluster