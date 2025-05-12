#!/bin/Rscript

# Script to convert accession numbers to taxonomy 
# Using the R package 'taxonomizr' (https://cran.r-project.org/web/packages/taxonomizr/readme/README.html)
library(taxonomizr)
library(tidyverse)

dir_data <- '/Users/lauragivens/Desktop/R/BZrookery_eDNA/Rdata'
dir_results <- "/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI/cutadapt/results"
setwd(dir_results)
#prepareDatabase('accessionTaxa.sql') #download data from NCBI and prepare SQLite database
## uses a LOT of hard drive space, FYI

######### Full database ######### 

blastResults<-read.table('dada2.uniques.BLAST.default.tsv',header=FALSE,stringsAsFactors=FALSE)
accessions<-strsplit(blastResults[,2],'\\|') #select the second column

taxaId<-accessionToTaxa(accessions,"~/accessionTaxa.sql")
taxResults <- getTaxonomy(taxaId,'~/accessionTaxa.sql')
taxResults <- as.data.frame(getTaxonomy(unique(taxaId),'~/accessionTaxa.sql'),getNames=TRUE) %>% #converts taxa ids to tax ranks
  rownames_to_column(.,"taxid") %>% mutate(.,taxid=str_trim(taxid))

t1d <- blastResults[1:2] #qseqid sseqid
t1d$taxid <- unlist(taxaId) %>% as.character()
t2d <- merge(t1d,taxResults,by.x='taxid',all.x=TRUE) %>% rename('qseqid'="V1","sseqid"="V2")


write.csv(taxaId,"taxonomizr.taxaID.csv")
write.csv(taxResults, "taxonomizr.taxResults.csv")

write.csv(t2d,'taxonomizr.merge.csv')
write_rds(t2d,'taxonomizr.merge.rds')



######### MAR database #########

marblastResults<-read.table('dada2.uniques.BLAST.martaxid.tsv',header=FALSE,stringsAsFactors=FALSE)
maraccessions<-strsplit(marblastResults[,2],'\\|') #pull accession numbers into a list

martaxaId<-accessionToTaxa(maraccessions,"~/accessionTaxa.sql") #converts accession number to tax id
martaxResults <- as.data.frame(getTaxonomy(unique(martaxaId),'~/accessionTaxa.sql')) %>% #converts taxa ids to tax ranks
                    rownames_to_column(.,"taxid") %>% mutate(.,taxid=str_trim(taxid))

t1 <- marblastResults[1:2] #qseqid sseqid
t1$taxid <- unlist(martaxaId) %>% as.character()
t2 <- merge(t1,martaxResults,by.x='taxid',all.x=TRUE) %>% rename('qseqid'="V1","sseqid"="V2")

write.csv(martaxaId,"taxonomizr.mar.taxaID.csv")
write.csv(martaxResults, "taxonomizr.mar.taxResults.csv")

write.csv(t2,'taxonomizr.mar.merge.csv')
write_rds(t2,'taxonomizr.mar.merge.rds')


save.image(paste0(dir_data,"/05_AccNo_Tax.RData"))
