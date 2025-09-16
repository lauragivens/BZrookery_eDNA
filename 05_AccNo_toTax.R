#!/bin/Rscript

########################### Setup ########################### 

# Script to convert accession numbers to taxonomy 
# Using the R package 'taxonomizr' (https://cran.r-project.org/web/packages/taxonomizr/readme/README.html)
library(taxonomizr)
library(tidyverse)

# set path variables  
dir_data <- '/Users/lauragivens/Desktop/R/BZrookery_eDNA/Rdata'
dir_results <- "/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI/cutadapt/results"
setwd(dir_results)
#prepareDatabase('accessionTaxa.sql') #download data from NCBI and prepare SQLite database
## uses a LOT of hard drive space, FYI

########################### Full database ########################### 

# read in BLAST results
# set header=FALSE because CLI doesn't output headers
# then make a vector of accession numbers  
# this will be used to get the NCBI taxid 
blastResults<-read.table('dada2.uniques.BLAST.wordsize.tsv',header=FALSE,stringsAsFactors=FALSE)
accessions<-strsplit(blastResults[,2],'\\|') #select the second column

# convert accessions to taxa 
# then get taxonomy from taxid 
# and merge taxids with tax ranks
taxaId<-accessionToTaxa(accessions,"~/accessionTaxa.sql")
taxResults <- getTaxonomy(taxaId,'~/accessionTaxa.sql')
taxResults <- as.data.frame(getTaxonomy(unique(taxaId),'~/accessionTaxa.sql'),getNames=TRUE) %>% 
  rownames_to_column(.,"taxid") %>% mutate(.,taxid=str_trim(taxid))

# select the qseqid and sseqid from BLAST output
# then merge the BLAST ids and taxaIDs into their own data frame
# this will result in a df (t2d) with 
# the full taxonomy (Kingdom-Species)
# the NCBI taxid 
# the dada2 sequence ID (qseqid)
# the NCBI hit sequence ID (sseqid)
t1d <- blastResults[1:2] #qseqid sseqid
t1d$taxid <- unlist(taxaId) %>% as.character()
t2d <- merge(t1d,taxResults,by.x='taxid',all.x=TRUE) %>% rename('qseqid'="V1","sseqid"="V2")

# save the taxid as a csv 
# save the taxonomy table with taxid as a csv 
# save the taxonomy table with metadata as a csv
write.csv(taxaId,"taxonomizr.wordsize.taxaID.csv")
write.csv(taxResults, "taxonomizr.wordsize.taxResults.csv")
write.csv(t2d,'taxonomizr.wordsize.merge.csv')

# save the taxonomy table with metadata as a rds object
write_rds(t2d,'taxonomizr.wordsize.merge.rds')

########################### MAR database ###########################

# read in BLAST results
# set header=FALSE because CLI doesn't output headers
# then make a vector of accession numbers  
# this will be used to get the NCBI taxid 
marblastResults<-read.table('dada2.uniques.BLAST.martaxid.tsv',header=FALSE,stringsAsFactors=FALSE)
maraccessions<-strsplit(marblastResults[,2],'\\|') #pull accession numbers into a list

# convert accessions to taxa 
# then get taxonomy from taxid 
# and merge taxids with tax ranks
martaxaId<-accessionToTaxa(maraccessions,"~/accessionTaxa.sql") 
martaxResults <- as.data.frame(getTaxonomy(unique(martaxaId),'~/accessionTaxa.sql')) %>% 
                    rownames_to_column(.,"taxid") %>% mutate(.,taxid=str_trim(taxid))

# select the qseqid and sseqid from BLAST output
# then merge the BLAST ids and taxaIDs into their own data frame
# this will result in a df (t2d) with 
# the full taxonomy (Kingdom-Species)
# the NCBI taxid 
# the dada2 sequence ID (qseqid)
# the NCBI hit sequence ID (sseqid)
t1 <- marblastResults[1:2] #qseqid sseqid
t1$taxid <- unlist(martaxaId) %>% as.character()
t2 <- merge(t1,martaxResults,by.x='taxid',all.x=TRUE) %>% rename('qseqid'="V1","sseqid"="V2")

# save the taxid as a csv 
# save the taxonomy table with taxid as a csv 
# save the taxonomy table with metadata as a csv
write.csv(martaxaId,"taxonomizr.mar.v2.taxaID.csv")
write.csv(martaxResults, "taxonomizr.mar.v2.taxResults.csv")
write.csv(t2,'taxonomizr.mar.merge.v2.csv')

# save the taxonomy table with metadata as a rds object
write_rds(t2,'taxonomizr.mar.merge.v2.rds')

########################### BOLD ###########################

# read in BOLD output  
# BOLD assigns 'no-match' to any taxonomy without a match
# replace with NA 
# then separate seqid into separate id and size columns  
bold_output <- readxl::read_xlsx(paste0(dir_results,"/dada2-uniqueseqs_identification_result.xlsx") )
bold_output[bold_output=='no-match'] <- NA
bold_output[bold_output=='IncompleteTaxonomy'] <- NA
bold_selected <- bold_output %>% mutate(.,seqid=word(id,sep=";"),.before=1) %>% as.data.frame()

# get a list of unique taxa identified to Species level  
# then get the taxid from those taxa
# and make a new dataframe associated species and taxid  
bold_spp <- bold_selected %>% filter(!is.na(Species)) %>% .$Species %>% unique()
bold_id <- getId(bold_spp,"~/accessionTaxa.sql",onlyScientific = TRUE)
bold_iddf <- cbind("Species"=bold_spp,"taxid"=bold_id)

t3 <- merge(bold_selected,bold_iddf,
            by="Species",
            all.x=TRUE,sort=FALSE) %>% 
  ungroup() %>% 
  relocate(c(Species,taxid),.after=Genus)

# save edited BOLD output as a new csv 
# save BOLD taxonomy with taxid as a new csv
# save BOLD taxonomy with taxid as an rds object
write.csv(bold_selected,paste0(dir_results,"/bold_identification.csv"))
write.csv(t3,paste0(dir_results,"/taxonomizr.bold.merge.rds"))
write_rds(t3,paste0(dir_results,"/taxonomizr.bold.merge.rds"))

save.image(paste0(dir_data,"/05_AccNo_Tax_wordsize.RData"))
