#Load packages
library(reticulate)
library(dada2)
library(data.table)
library(tidyverse)
library(Biostrings)

# activate conda environment
use_condaenv(condaenv = 'cutadapt', required=TRUE)

# set directories
setwd('/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI')
dir_home <- getwd()
dir_raw<-paste0(dir_home,'/fastq')
dir_cut<-paste0(dir_home,'/cutadapt')
if(!dir.exists(dir_cut)) dir.create(dir_cut)
  
cutadapt<-"/Users/lauragivens/miniconda3/envs/cutadapt/bin/cutadapt"
# primer sequences
# sequences and citations backtracked from SERC metabarcoding digital notebook 
# citation:https://frontiersinzoology.biomedcentral.com/articles/10.1186/1742-9994-10-34
FwdPrimer=c("GGWACWGGWTGAACWGTWTAYCCYCC") #ILL-mlLCOF1
FwdPrimerRC=dada2:::rc(FwdPrimer) 
#RevPrimer=c("TAIACYTCIGGRTGICCRAARAAYCA") #ILL-igHCO2198R #rc does not recognzie I
RevPrimer=c("TANACYTCNGGRTGNCCRAARAAYCA") #ILL-jgHCO2198R
RevPrimerRC=dada2:::rc(RevPrimer)

# Step 1: Demux: completed by LAB

# Step 2: Remove primers  
# get sample names from files (assuming output follows naming convention of NAME_R1_001.fastq.gz, NAME_R1_001.fastq.gz)
fnFs <- sort(list.files(dir_raw, pattern="R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(dir_raw, pattern="R2_001.fastq.gz", full.names = TRUE))

# make first subdir for trimmed samples
dir_trim <- file.path(dir_cut,'trim')
if(!dir.exists(dir_trim)) dir.create(dir_trim)

# list out file names for each sample to go into cutadapt 
fnFs_trim <- file.path(dir_trim,basename(fnFs))
fnRs_trim <- file.path(dir_trim,basename(fnRs))

# combine F/R primers and their RC for R1 and R2 (forward and reverse reads)
R1.flags <- paste('-g',FwdPrimer,
                  '-a',RevPrimerRC)
R2.flags <- paste('-G',FwdPrimerRC,
                  '-A',RevPrimer)

# run cutadapt for each sample and write to dir_trim
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, 
                             R2.flags,
                             "--discard-untrimmed",
                             "--json",paste0(fnFs[i],".cutadapt.json"),
                             "--max-n 0", #discards reads containing more than  0 N (ambiguous) bases
                             #--overlap 3, #may need to change the minimum overlap length between read and the adapter sequence. default is 3 bases
                             "-o", fnFs_trim[i], #out.1.fastq
                             "-p", fnRs_trim[i], #out.2.fastq 
                             fnFs[i], #input forward files
                             fnRs[i] #input reverse files
                             ))
}


system2(cutadapt, args = c(R1.flags, 
                           R2.flags,
                           "--discard-untrimmed",
                           "--verbose",
                           "--json",paste0(fnFs[1],".cutadapt.json"),
                           "--max-n 0", #discards reads containing more than  0 N (ambiguous) bases
                           #--overlap 3, #may need to change the minimum overlap length between read and the adapter sequence. default is 3 bases
                           "-o", fnFs_trim[1], #out.1.fastq
                           "-p", fnRs_trim[1], #out.2.fastq 
                           fnFs[1], #input forward files
                           fnRs[1] #input reverse files
))

Sys.Date()
sessionInfo()
