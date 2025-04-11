library(Biostrings)
library(lulu)

# cluster into ASVs
## https://github.com/tobiasgf/lulu

# requires a table
# samples as columns, seqs as rows
# unique seq id as rownames

# requires a fasta file of sequences

###########################################

dir_results <- "/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI/cutadapt/results"
uniqueseqs <- readDNAStringSet(paste0(dir_results,'/dada2-uniqueseqs.fasta'))
seqtab.nochim <- readRDS(file=paste0(dir_results,'/seqtab.nochim.rds'))
                         
# convert seqtab
lulu.seqtab.nochim <- as.data.frame(seqtab.nochim)
colnames(lulu.seqtab.nochim) <- word(names(uniqueseqs),sep=';')
rownames(lulu.seqtab.nochim) <- word(rownames(lulu.seqtab.nochim),1,4,sep='_')
lulu.seqtab.nochim <- t(lulu.seqtab.nochim) %>% as.data.frame()
#lulu.seqtab.nochim$OTUID <- rownames(lulu.seqtab.nochim)
lulu.seqtab.nochim[1:5,1:5]

# convert fasta file
lulu.uniqueseqs <- uniqueseqs
names(lulu.uniqueseqs) <- word(names(lulu.uniqueseqs),sep=';')
head(names(lulu.uniqueseqs),2)

# 2. produce a match list
## essentially pairwise matching of sequences

### doing this outside of R
## First produce a blastdatabase with the OTUs
# makeblastdb -in dada2-uniqueseqs.fasta -parse_seqids -dbtype nucl
## Then blast the OTUs against the database
# blastn -db dada2-uniqueseqs.fasta -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query dada2-uniqueseqs.fasta

## read in 
match_list <- read_delim(paste0(dir_results,'/match_list.txt'),col_names=F,delim='\t')
match_list$X1 <- word(match_list$X1,sep=';')
match_list$X2 <- word(match_list$X2,sep=';')
match_list[1:3,1:3]

# 3. Run LULU curation 
curated_result <- lulu(lulu.seqtab.nochim, match_list)

# Save  
saveRDS(curated_result,file=paste0(dir_results,"/lulu-clustertable.rds"))
save.image("/Users/lauragivens/Desktop/R/BZrookery_eDNA/Rdata/03_LULU.RData")
