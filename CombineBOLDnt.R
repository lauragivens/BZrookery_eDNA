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

taxa_nt_sub <- readRDS(paste0(dir_results,"/taxtable.v2.nt.rds"))
taxa_nt_troph <- readRDS(paste0(dir_results,"/taxtable.v2.nt.wtroph.rds"))

bold_sub <- bold_selected %>% filter(!is.na(Phylum) & Phylum!="IncompleteTaxonomy")

View(bold_sub)

combined_tax <- rbind((bold_sub[c(1,3:9,11)]  %>% 
                         mutate(db=1)),
                      (taxa_nt_sub[2:7] %>% 
                         rownames_to_column(var="seqid") %>%
                         mutate(db=2,pct_identity=NA,records=NA))
                      ) 

# select duplicates
# combined_tax %>% 
#  group_by(seqid) %>% filter(n() == 2) %>%
#  View()

check.tax <- combined_tax %>% 
  group_by(seqid) %>% filter(n() == 2) %>% 
  arrange(seqid) %>%
  # mutate(across(everything(),!duplicated())) %>% 
  mutate(across(where(is.character), ~ duplicated(.x), .names = "{.col}_dup"),
         across(where(is.logical), ~ case_when( (.x==TRUE) ~ .x))
  ) 

check.tax %>% 
  group_by(seqid) %>%
  summarize(across(where(is.character),~ unique(.x))) %>%
  View()

combined_tax %>% 
  group_by(seqid) %>% 
  filter(n()==2) %>% 
  summarise(across(where(is.character), ~ na.omit(.x))) %>% View()

