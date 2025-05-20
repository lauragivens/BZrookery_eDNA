library(phyloseq)
library(vegan)
library(tidyverse)

dir_data <- '/Users/lauragivens/Desktop/R/BZrookery_eDNA/Rdata'
dir_man <- "/Volumes/Fuji/Mangroves"
dir_results <- "/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI/cutadapt/results"
setwd(dir_results)

bold_output <- readxl::read_xlsx("dada2-uniqueseqs_identification_result.xlsx") 
bold_output[bold_output=='no-match'] <- NA
bold_selected <- bold_output %>% mutate(.,seqid=word(id,sep=";"),.before=1) %>% as.data.frame()
rownames(bold_selected) <- bold_selected$seqid

# upload asv table  
curated_lulu <- readRDS('lulu-clustertable.rds')
curated_asv <- curated_lulu$curated_table
cols_asv <- word(colnames(curated_asv),sep = "_",end=3) 
colnames(curated_asv) <- cols_asv

# upload metadata
samplelist <- read.csv("/Users/lauragivens/Desktop/R/BZrookery_eDNA/Metadata_BZrookery.csv")
rownames(samplelist) <- samplelist$SampleName_Long
samplelist$SampleName_Long <- NULL

# assemble ps object
otu <- otu_table((curated_asv),taxa_are_rows = TRUE)
meta <- sample_data(samplelist)
tax <- tax_table(as.matrix(bold_selected))

ps <- phyloseq(otu,meta,tax)
ps

write.csv(bold_selected,paste0(dir_results,"/taxtable.bold.csv"))
saveRDS(bold_selected,paste0(dir_results,"/taxtable.bold.rds"))
saveRDS(ps,paste0(dir_results,"/ps.bold.rds"))


save.image(paste0(dir_data,'/06_BOLDtoPhyloseq.RData'))
