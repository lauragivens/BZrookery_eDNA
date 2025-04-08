#!/bin/Rscript

# Script to convert accession numbers to taxonomy 
# Using the R package 'taxonomizr' (https://cran.r-project.org/web/packages/taxonomizr/readme/README.html)
library(taxonomizr)

prepareDatabase('accessionTaxa.sql') #download data from NCBI and prepare SQLite database
## uses a LOT of hard drive space, FYI

blastResults<-read.table('dada2.uniques.BLAST.default.tsv',header=FALSE,stringsAsFactors=FALSE)
accessions<-strsplit(blastResults[,2],'\\|') #select the second column

taxaId<-accessionToTaxa(accessions,"accessionTaxa.sql")
taxResults <- getTaxonomy(taxaId,'accessionTaxa.sql')

write.csv(taxaId,"taxonomizr.taxaID.csv")
write.csv(taxResults, "taxonomizr.taxResults.csv")


#########

marblastResults<-read.table('dada2.uniques.BLAST.martaxid.tsv',header=FALSE,stringsAsFactors=FALSE)
maraccessions<-strsplit(marblastResults[,2],'\\|')

martaxaId<-accessionToTaxa(maraccessions,"accessionTaxa.sql")
martaxResults <- getTaxonomy(martaxaId,'accessionTaxa.sql')

write.csv(martaxaId,"taxonomizr.mar.taxaID.csv")
write.csv(martaxResults, "taxonomizr.mar.taxResults.csv")
