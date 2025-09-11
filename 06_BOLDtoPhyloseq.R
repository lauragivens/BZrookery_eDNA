#!/bin/Rscript

########################### Setup ########################### 

library(phyloseq)
library(vegan)
library(tidyverse)

# set path variables  
dir_data <- '/Users/lauragivens/Desktop/R/BZrookery_eDNA/Rdata'
dir_man <- "/Volumes/Fuji/Mangroves"
dir_results <- "/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI/cutadapt/results"
setwd(dir_results)

###########################  Read in files ########################### 

# upload BOLD taxonomy results 
# Use BOLD taxonomy with NCBI taxids  
# then assign seqid as rownames for phyloseq  
bold_tax <- read_rds(paste0(dir_results,"/taxonomizr.bold.merge.rds"))
#bold_selected <- read.csv(paste0(dir_results,'/bold_identification.csv')) %>% select(-X)
rownames(bold_tax) <- bold_tax$seqid
bold_tax <- bold_tax %>% select(-c(seqid,id))

# upload asv table from LULU    
# select only the curated table  
# then make a vector of the column names and truncate to match metadata  
# assign that vector as the new column names
curated_lulu <- readRDS('lulu-clustertable.rds')
curated_asv <- curated_lulu$curated_table
cols_asv <- word(colnames(curated_asv),sep = "_",end=3) 
colnames(curated_asv) <- cols_asv

# upload metadata file  
# assign the sample names as rownames 
# then remove from metadata
samplelist <- read.csv("/Users/lauragivens/Desktop/R/BZrookery_eDNA/Metadata_BZrookery.csv")
rownames(samplelist) <- samplelist$SampleName_Long
samplelist$SampleName_Long <- NULL

########################### Assemble ps object ########################### 
otu <- otu_table((curated_asv),taxa_are_rows = TRUE)
meta <- sample_data(samplelist)
tax <- tax_table(as.matrix(bold_tax))

ps <- phyloseq(otu,meta,tax)
ps

########################### Save ########################### 
write.csv(bold_tax,paste0(dir_results,"/taxtable.bold.taxid.csv"))
saveRDS(bold_tax,paste0(dir_results,"/taxtable.bold.taxid.rds"))
saveRDS(ps,paste0(dir_results,"/ps.bold.taxid.rds"))


save.image(paste0(dir_data,'/06_BOLDtoPhyloseq.v2.RData'))
