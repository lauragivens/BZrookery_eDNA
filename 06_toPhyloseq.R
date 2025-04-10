library(phyloseq)
library(vegan)
library(tidyverse)

dir <- "/Volumes/Fuji/Mangroves"
dir_results <- "/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI/cutadapt/results"
setwd(dir_results)

# upload BLAST files 
blast_output <- read.delim('dada2.uniques.BLAST.martaxid.tsv',header = FALSE)
colnames(blast_output) <- c('qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore')

# upload accession - to - taxid files
taxid_mar <- read.csv('taxonomizr.mar.taxaID.csv')
taxresults_mar <- read.csv('taxonomizr.mar.taxResults.csv')
t2 <- readRDS('taxonomizr.mar.merge.rds')

# combine taxa names and blast output 
taxa_df <- left_join(blast_output,t2,by=c('qseqid','sseqid'))

# choose first occurrence of each asv 
#qseqid_unique <-taxa_df$qseqid[!duplicated(taxa_df$qseqid)]
qseqid_unique <- taxa_df[!duplicated(taxa_df['qseqid']),]

# use %age match to truncate assignments at appropriate rank
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

# upload trophic data
bzfishmeta <- read.csv(paste0(dir,'/BelizeanFishSpecies.csv'),header=TRUE)
bzfishmeta <- bzfishmeta %>% rename(., 'taxid'="NCBI_taxid") %>% mutate(taxid=as.character(taxid)) %>% select(-c(starts_with("AccNo"),'Species.name'))
bzfishmeta <- bzfishmeta %>% filter(!is.na(taxid))

# combine trophic and tax data  
taxa_troph <- left_join((taxa_final %>% select(c(seqid,14:21))),bzfishmeta,by='taxid')
rownames(taxa_troph) <- taxa_troph$seqid
taxa_troph <- taxa_troph %>% select(-c('seqid')) %>% relocate(taxid,.after='species')

# upload asv table  
curated_lulu <- readRDS('lulu-clustertable.rds')
curated_asv <- curated_lulu$curated_table
cols_asv <- word(colnames(curated_asv),sep = "_",end=3) 
colnames(curated_asv) <- cols_asv

# upload metadata
samplelist <- read.csv(paste0(dir,"/TotalSampleList-BZ.csv"))
rownames(samplelist) <- samplelist$X
samplelist$X <- NULL

# assemble ps object
otu <- otu_table((curated_asv),taxa_are_rows = TRUE)
meta <- sample_data(samplelist)
tax <- tax_table(as.matrix(taxa_sub))
taxtroph <- tax_table(as.matrix(taxa_troph))

ps <- phyloseq(otu,meta,tax)
ps.troph <- phyloseq(otu,meta,taxtroph)
ps
ps.troph

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
