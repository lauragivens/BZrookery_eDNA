library(phyloseq)
library(vegan)
library(tidyverse)

dir_data <- '/Users/lauragivens/Desktop/R/BZrookery_eDNA/Rdata'
dir_home <- '/Users/lauragivens/Downloads/Demo'
dir_results <- "/Users/lauragivens/Downloads/Demo/cutadapt/results"
setwd(dir_results)

# upload BLAST files 
blast_output <- read.delim('dada2.uniques.BLAST.default.tsv',header = FALSE)
colnames(blast_output) <- c('qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore')

# upload accession - to - taxid files
taxid <- read.csv('taxonomizr.v2.taxaID.csv')
taxresults <- read.csv('taxonomizr.v2.taxResults.csv')
t2d <- readRDS('taxonomizr.v2.merge.rds')
names(t2d)[4:10] <- str_to_title(names(t2d)[4:10])
t2dd <- distinct(t2d)

# combine taxa names and blast output 
taxa_df <- left_join(blast_output,t2dd,by=c('qseqid','sseqid')) 

# choose first occurrence of each asv 
qseqid_unique <- taxa_df[!duplicated(taxa_df['qseqid']),]

##########################################################################################################
# use %age match to truncate assignments at appropriate rank (blast nt database)
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
taxa_final <- taxa_final %>% mutate(Species=case_when(str_detect(Species,"sp.") ~ NA,
                                                      str_detect(Species,"spp.") ~ NA,
                                                      str_detect(Species,"cf.") ~ NA,
                                                      str_detect(Species,"uncultur") ~ NA,
                                                      str_detect(Species,"environ") ~ NA,
                                                      str_detect(Species,"acterium") ~ NA,
                                                      str_count(Species," ") != 1 ~ NA,
                                                      .default=Species)
) 

rownames(taxa_final) <- taxa_final$seqid
taxa_sub <- select(taxa_final, c(15:21))


##########################################################################################################
# upload asv table  
curated_lulu <- readRDS('lulu-clustertable.rds')
curated_asv <- curated_lulu$curated_table
cols_asv <- word(colnames(curated_asv),sep = "_S",1) 
colnames(curated_asv) <- cols_asv

# upload metadata
samplelist <- read.csv(paste0(dir_home,"/Metadata_BZrookery.csv"))
rownames(samplelist) <- samplelist$SampleName_Long
samplelist$SampleName_Long <- NULL

meta.tab <- samplelist[rownames(samplelist)%in%cols_asv,]

# assemble ps object
otu <- otu_table((curated_asv),taxa_are_rows = TRUE)
meta <- sample_data(meta.tab)
tax <- tax_table(as.matrix(taxa_sub))

ps <- phyloseq(otu,meta,tax)
ps

##########################################################################################################
# Save nt 
saveRDS(ps,paste0(dir_results,"/ps.rds"))

saveRDS(taxa_sub,paste0(dir_results,"/taxtable.rds"))
write.csv(taxa_sub,paste0(dir_results,"/taxtable.csv"))

save.image(paste0(dir_data,'/06_toPhyloseq.RData'))
