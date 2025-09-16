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

# read in BOLD file 
# seqid as rownames
# has NCBI taxid included
bold_output <- readRDS(paste0(dir_results,"/taxtable.bold.taxid.rds"))
bold_sub <- bold_output %>% filter(!is.na(Phylum))

# read in taxonomy tables 
taxa_nt_sub <- readRDS(paste0(dir_results,"/taxtable.wordsize.nt.rds"))
taxa_nt_troph <- readRDS(paste0(dir_results,"/taxtable.wordsize.nt.wtroph.rds"))

###########################  Merge NCBI and BOLD results ########################### 

# make vectors of the column names 
# for bold and ncbi output  
bold.cols <- colnames(bold_sub)
nt.cols <- colnames(taxa_nt_sub)

# make a bold dataframe 
# with a new column for rownames
# and fill out the df with additional columns that were in the ncbi df  
bold.tomerge <- (bold_sub %>% 
                   rownames_to_column(var="seqid") %>% 
                   mutate(db=1)) 
bold.tomerge[,nt.cols[!nt.cols%in%bold.cols]] <- NA

# make a ncbi dataframe 
# with a new column for rownames
# and fill out the df with additional columns that were in the bold df  
nt.tomerge <- (taxa_nt_sub %>% 
                 rownames_to_column(var="seqid") %>%
                 mutate(db=2))
nt.tomerge[,bold.cols[!bold.cols%in%nt.cols]] <- NA

# merge the dfs 
# add an additional column that identifies whether a tax assignment is from BOLD or NCBI  
combined_tax <- rbind(bold.tomerge,
                      nt.tomerge
                      ) %>% mutate(db=case_when(db==1 ~ "BOLD",
                                                db==2 ~ "nt")) %>% 
  relocate(Domain,.before=Phylum) %>% 
  arrange(seqid)

########################### Compare tax assignments and truncate ########################### 

# we are going to take the combined dataframe 
# which has one row per assignment per database  
# which means that some seqids may have two rows 
# we need to get to one row per seqid  
tax.combined.df <- combined_tax %>%   
  group_by(seqid) %>% 
  # for each sequence id, we are going to summarize every column with summarize_all
  # ignoring any NA values  
  # we keep the unique values  
  # toString() essentially works like paste() 
  # it collapses everything in a column into a single string separated by ","
  # so for a seqid, if the tax assignment for a level is the same, it will output the unique value (e.g., Chordata)
  # but if BOLD and NCBI disagree, it will output a comma-separated value (e.g., Chordata,Mollusca)
  summarise_all(list(~toString(unique(na.omit(.))))) %>% #unique for cells that are the same / duplicates
  # we replace all the blank values with NA 
  # so we can use the helpful tidy functions with NA
  mutate_all(., list(~na_if(.,""))) %>%
  # we then check the higher level taxa for commas    
  # if we detect a comma in a higher level taxa
  # but no comma in Species
  # we are going to assume that that is a disagreement in naming convention 
  # (one db uses an older/unaccepted name or assigned taxa to a sublevel)
  # and we select the first value  
  mutate(across(c(Phylum,Class,Order,Family,Genus), #check the higher level taxa
                ~ case_when(str_detect(.,",") & !str_detect(Species,",") & !is.na(Species) ~ word(.,sep=","), 
                            .default=.x)
                )
         ) %>% 
  # now we make some decisions about where to truncate assignments 
  # based on disagreements between the databases  
  # (indicated by a comma)
  # basic structure is 
  # for each level, if there is a comma and the next lower level is NA
  # that level is assigned to NA as well 
  # but if the next lower level is NOT NA  
  # truncate the current level to the first assignment  
  # 
  # this essentially is assuming that a comma at higher levels surrounded by agreement at lower levels 
  # indicates a difference in naming convention 
  # instead of disagreement between actual assignments  
  mutate(Species=case_when(str_detect(Species,",") ~ NA,
                           #if there is a comma in species, we replace that with NA
                           .default=Species),
         Genus=case_when(str_detect(Genus,",") & is.na(Species) ~ NA, 
                         #if Genus has a comma and Species is NA 
                         #Genus is NA
                         str_detect(Genus,",") & !is.na(Species) ~ word(Genus,sep=","),
                         #if Genus has a comma and Species has a value
                         #keep the first value of Genus
                         .default=Genus),
         Family=case_when(str_detect(Family,",") & is.na(Genus) ~ NA, 
                          #if Family has a comma and Genus is NA 
                          #Family is NA
                          str_detect(Family,",") & !is.na(Genus) ~ word(Family,sep=","),
                          #if Family has a comma and Genus has a value
                          #keep the first value of Family
                          .default=Family),
         Order=case_when(str_detect(Order,",") & is.na(Family) ~ NA, 
                         #if Order has a comma and Family is NA 
                         #Order is NA
                         str_detect(Order,",") & !is.na(Family) ~ word(Order,sep=","),
                         #if Order has a comma and Family has a value
                         #keep the first value of Order
                         .default=Order),
         Class=case_when(str_detect(Class,",") & is.na(Order) ~ NA, 
                         #if Class has a comma and Order is NA 
                         #Class is NA
                         str_detect(Class,",") & !is.na(Order) ~ word(Class,sep=","),
                         #if Class has a comma and Order has a value
                         #keep the first value of Class
                         .default=Class),
         Phylum=case_when(str_detect(Phylum,",") & is.na(Class) ~ NA, 
                          #if Phylum has a comma and Class is NA 
                          #Phylum is NA
                          str_detect(Phylum,",") & !is.na(Class) ~ word(Phylum,sep=","),
                          #if Phylum has a comma and Class has a value
                          #keep the first value of Phylum
                          .default=Phylum)
  ) 

#how many unique Chordate species  
tax.combined.df %>% filter(Phylum=="Chordata") %>% .$Species %>% unique() %>% length()
#are any seqids repeated? 
summary(duplicated(tax.combined.df$seqid))

########################### Replace middle NA values ########################### 

# Now we have some weird values from the nt tax database that didn't fill out correctly
# they leave an NA and fill the rest of the taxonomy 
# how do we fill that in if we have the information  

# TBD  


########################### Add trait info ########################### 

# upload trophic data
bzfishmeta <- read.csv(paste0(dir_man,'/BelizeanFishSpecies.csv'),header=TRUE)
bzfishmeta <- bzfishmeta %>% dplyr::rename(., "taxid" = "NCBI_taxid", 
                                           "Species.Abundance" = "Abundance") %>% 
  mutate(taxid=as.character(taxid)) %>% 
  select(-c(starts_with("AccNo"),'Species.name')) 

names(bzfishmeta) <- trimws(names(bzfishmeta)) 
names(bzfishmeta) <- gsub("_",".",names(bzfishmeta))

bzfishmeta <- bzfishmeta %>% filter(!is.na(taxid))
bzfishmeta <- bzfishmeta %>% mutate(Trophic.Level=word(Trophic.Index,sep=" ") %>% as.numeric(),
                                    Trophic.Rounded=round(Trophic.Level,digits=0),
                                    .after="Trophic.Index"
)

# combine trophic and tax data  
taxa_troph <-  left_join((tax.combined.df #%>% select(c(seqid,14:21))
             ),
            bzfishmeta,by='taxid') %>% as.data.frame()

rownames(taxa_troph) <- taxa_troph$seqid
taxa_troph <- taxa_troph %>% select(-c('seqid')) %>% relocate(taxid,.after='Species')

########################### Make new PS object ########################### 
# conform tax table  
tax.combined.df <- as.data.frame(tax.combined.df)
rownames(tax.combined.df) <- tax.combined.df$seqid

# upload asv table  
curated_lulu <- readRDS(paste0(dir_results,'/lulu-clustertable.rds'))
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
tax <- tax_table(as.matrix(tax.combined.df %>% select(-seqid)))
tax.troph <- tax_table(as.matrix(taxa_troph))

ps <- phyloseq(otu,meta,tax)
ps.troph <- phyloseq(otu,meta,tax.troph)
ps
ps.troph

########################### Save ########################### 
saveRDS(ps,paste0(dir_results,"/ps.combined.wordsize.rds"))
saveRDS(ps.troph,paste0(dir_results,"/ps.troph.combined.wordsize.rds"))

write.csv(tax.combined.df,paste0(dir_results,"/tax.combined.wordsize.df.csv"))
saveRDS(tax.combined.df,paste0(dir_results,"/tax.combined.wordsize.df.rds"))

write.csv(taxa_troph,paste0(dir_results,"/tax.troph.combined.wordsize.df.csv"))
saveRDS(taxa_troph,paste0(dir_results,"/tax.troph.combined.wordsize.df.rds"))

save.image(paste0(dir_data,'/CombineBOLDnt.wordsize.RData'))





  
