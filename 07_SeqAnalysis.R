library(phyloseq)
library(tidyverse)

dir_data <- '/Users/lauragivens/Desktop/R/BZrookery_eDNA/Rdata'
dir_results <- "/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI/cutadapt/results"
setwd(dir_results)

# Load 
ps <- readRDS(paste0(dir_results,"/ps.rds"))
ps.troph <- readRDS(paste0(dir_results,"/ps.troph.rds"))

taxa <- readRDS(paste0(dir_results,"/taxtable.rds"))
taxa_troph <- readRDS(paste0(dir_results,"/taxtable.wtroph.rds"))
curated_lulu <- readRDS('lulu-clustertable.rds')
curated_asv <- readRDS(paste0(dir_results,"/asvtable.rds"))
samplelist <- readRDS(paste0(dir_results,'/metadata.rds'))

