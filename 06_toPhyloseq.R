library(phyloseq)
library(vegan)
library(tidyverse)

dir_data <- '/Users/lauragivens/Desktop/R/BZrookery_eDNA/Rdata'
dir_man <- "/Volumes/Fuji/Mangroves"
dir_results <- "/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI/cutadapt/results"
setwd(dir_results)

# upload BLAST files 
blast_output <- read.delim('dada2.uniques.BLAST.martaxid.tsv',header = FALSE)
colnames(blast_output) <- c('qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore')

blast_nt_output <- read.delim('dada2.uniques.BLAST.default.tsv',header = FALSE)
colnames(blast_nt_output) <- c('qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore')

# upload accession - to - taxid files
taxid_mar <- read.csv('taxonomizr.mar.taxaID.csv')
taxresults_mar <- read.csv('taxonomizr.mar.taxResults.csv')
t2 <- readRDS('taxonomizr.mar.merge.rds')
names(t2)[4:10] <- str_to_title(names(t2)[4:10])

taxid <- read.csv('taxonomizr.taxaID.csv')
taxresults <- read.csv('taxonomizr.taxResults.csv')
t2d <- readRDS('taxonomizr.merge.rds')
names(t2d)[4:10] <- str_to_title(names(t2d)[4:10])
t2dd <- distinct(t2d)

# combine taxa names and blast output 
taxa_df <- left_join(blast_output,t2,by=c('qseqid','sseqid'))
taxa_nt_df <- left_join(blast_nt_output,t2dd,by=c('qseqid','sseqid')) #gets a warning if I don't manipulate t2d; #warning: Detected an unexpected many-to-many relationship between `x` and `y`.

# choose first occurrence of each asv 
#qseqid_unique <-taxa_df$qseqid[!duplicated(taxa_df$qseqid)]
qseqid_unique <- taxa_df[!duplicated(taxa_df['qseqid']),]
qseqid_nt_unique <- taxa_nt_df[!duplicated(taxa_nt_df['qseqid']),]

##########################################################################################################
# use %age match to truncate assignments at appropriate rank (mar fish taxa)
for (i in 1:dim(qseqid_unique)[1]) {
  seqs <- qseqid_unique[i,][[1]] #get the unique seq id 
  check <- taxa_df[taxa_df$qseqid==seqs,] #subset taxa_df to include only those seqids 
  truncated <- check[1,] #add the top hit as a separate variable
  rows <- nrow(check) #get the number of rows 
  
  diffs <- abs(check[1,]$pident-check[rows,]$pident) #get the absolute value of differences in pident between top and last hit

    if (diffs > 0) { #if there is a difference, do
      tona <- ifelse((diffs>0 & diffs<5), 20, #truncate to genus if pident>95
                     if_else((diffs>5 & diffs <13), 19:20, #family
                             if_else((diffs>13 & diffs < 17), 18:20, #order
                                     if_else((diffs>17 & diffs <19), 17:20, #class
                                             if_else((diffs>19 & diffs <21), 16:20, #phylum
                                                     if_else((diffs>21 & diffs <29), 5:20, #kingdom
                                                             14:20 #remove all classifications
                                                     ))))))
    
    truncated[1,][tona] <- NA #replace the given index with NA
    
  }
  else { #otherwise, keep the first row
    truncated <- truncated #could probably also tell it qseqid_unique[i,] <- qseqid_unique[i,]
    # check that 
  }
  
  qseqid_unique[i,] <- truncated
}
taxa_final <- qseqid_unique %>% mutate(.,seqid=word(qseqid,sep=";"),.before=1) 
rownames(taxa_final) <- taxa_final$seqid
taxa_sub <- select(taxa_final, c(15:21))

##########################################################################################################
# use %age match to truncate assignments at appropriate rank (blast nt database)
for (i in 1:dim(qseqid_nt_unique)[1]) {
  seqs_nt <- qseqid_nt_unique[i,][[1]] #get the unique seq id 
  check_nt <- taxa_nt_df[taxa_nt_df$qseqid==seqs_nt,] #subset taxa_df to include only those seqids 
  truncated_nt <- check_nt[1,] #add the top hit as a separate variable
  rows_nt <- nrow(check_nt) #get the number of rows 
  
  diffs_nt <- abs(check_nt[1,]$pident-check_nt[rows_nt,]$pident) #get the absolute value of differences in pident between top and last hit
  
  if (diffs_nt > 0) { #if there is a difference, do
    tona <- ifelse((diffs_nt>0 & diffs_nt<5), 20, #truncate to genus if pident>95
                   if_else((diffs_nt>5 & diffs_nt <13), 19:20, #family
                           if_else((diffs_nt>13 & diffs_nt < 17), 18:20, #order
                                   if_else((diffs_nt>17 & diffs_nt <19), 17:20, #class
                                           if_else((diffs_nt>19 & diffs_nt <21), 16:20, #phylum
                                                   if_else((diffs_nt>21 & diffs_nt <29), 5:20, #kingdom
                                                           14:20 #remove all classifications
                                                   ))))))
    
    truncated_nt[1,][tona] <- NA #replace the given index with NA
    
  }
  else { #otherwise, keep the first row
    truncated_nt <- truncated_nt #could probably also tell it qseqid_nt_unique[i,] <- qseqid_nt_unique[i,]
    # check that 
  }
  
  qseqid_nt_unique[i,] <- truncated_nt
}
taxa_nt_final <- qseqid_nt_unique %>% mutate(.,seqid=word(qseqid,sep=";"),.before=1) 
rownames(taxa_nt_final) <- taxa_nt_final$seqid
taxa_nt_sub <- select(taxa_nt_final, c(15:21))

##########################################################################################################
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
taxa_troph <- left_join((taxa_final %>% select(c(seqid,14:21))),bzfishmeta,by='taxid')
rownames(taxa_troph) <- taxa_troph$seqid
taxa_troph <- taxa_troph %>% select(-c('seqid')) %>% relocate(taxid,.after='Species')

taxa_nt_troph <- left_join((taxa_nt_final %>% select(c(seqid,14:21))),bzfishmeta,by='taxid')
rownames(taxa_nt_troph) <- taxa_nt_troph$seqid
taxa_nt_troph <- taxa_nt_troph %>% select(-c('seqid')) %>% relocate(taxid,.after='Species')

##########################################################################################################
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
tax <- tax_table(as.matrix(taxa_sub))
taxtroph <- tax_table(as.matrix(taxa_troph))
taxnt <- tax_table(as.matrix(taxa_nt_sub))
taxnttroph <-tax_table(as.matrix(taxa_nt_troph))

ps <- phyloseq(otu,meta,tax)
ps.troph <- phyloseq(otu,meta,taxtroph)
ps.nt <- phyloseq(otu,meta,taxnt)
ps.nt.troph <- phyloseq(otu,meta,taxnttroph)
ps
ps.troph
ps.nt
ps.nt.troph

##########################################################################################################
# Save mar
saveRDS(ps,paste0(dir_results,"/ps.rds"))
saveRDS(ps.troph,paste0(dir_results,"/ps.troph.rds"))

saveRDS(taxa_sub,paste0(dir_results,"/taxtable.rds"))
saveRDS(taxa_troph,paste0(dir_results,"/taxtable.wtroph.rds"))
saveRDS(curated_asv,paste0(dir_results,"/asvtable.rds"))
saveRDS(samplelist,paste0(dir_results,'/metadata.rds'))


write.csv(taxa_sub,paste0(dir_results,"/taxtable.csv"))
write.csv(taxa_troph,paste0(dir_results,"/taxtable.wtroph.csv"))
write.csv(curated_asv,paste0(dir_results,"/asvtable.csv"))
write.csv(samplelist,paste0(dir_results,'/metadata.csv'))

# Save nt 
saveRDS(ps.nt,paste0(dir_results,"/ps.nt.rds"))
saveRDS(ps.nt.troph,paste0(dir_results,"/ps.nt.troph.rds"))

saveRDS(taxa_nt_sub,paste0(dir_results,"/taxtable.nt.rds"))
saveRDS(taxa_nt_troph,paste0(dir_results,"/taxtable.nt.wtroph.rds"))

write.csv(taxa_nt_sub,paste0(dir_results,"/taxtable.nt.csv"))
write.csv(taxa_nt_troph,paste0(dir_results,"/taxtable.nt.wtroph.csv"))




save.image(paste0(dir_data,'/06_toPhyloseq.RData'))
