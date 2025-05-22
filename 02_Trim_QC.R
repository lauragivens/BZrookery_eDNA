#Load packages
library(reticulate)
library(dada2)
library(data.table)
library(tidyverse)
library(Biostrings)

# activate conda environment
use_condaenv(condaenv = 'cutadapt', required=TRUE)

# set directories
setwd('/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI')
dir_home <- getwd() #current working directory  
dir_raw<-paste0(dir_home,'/fastq') #folder where raw fastq files are located
dir_cut<-paste0(dir_home,'/cutadapt') #folder to write cutadapt results to 
if(!dir.exists(dir_cut)) dir.create(dir_cut)  #create the cutadapt directory if it does not exist  

dir_data <- '/Users/lauragivens/Desktop/R/BZrookery_eDNA/Rdata' #folder to save Rdata files to 
if(!dir.exists(dir_data)) dir.create(dir_data) #create the folder if it does not exist 
  
cutadapt<-"/Users/lauragivens/miniconda3/envs/cutadapt/bin/cutadapt" #location of cutadapt on the machine  
# primer sequences
# sequences and citations backtracked from SERC metabarcoding digital notebook 
# citation:https://frontiersinzoology.biomedcentral.com/articles/10.1186/1742-9994-10-34
FwdPrimer=c("GGWACWGGWTGAACWGTWTAYCCYCC") #ILL-mlLCOF1
FwdPrimerRC=dada2:::rc(FwdPrimer) 

RevPrimer=c("TANACYTCNGGRTGNCCRAARAAYCA") #ILL-jgHCO2198R
RevPrimerRC=dada2:::rc(RevPrimer)

# Step 1: Demux: completed by LAB

# Step 2: Remove primers  
# get sample names from files (assuming output follows naming convention of NAME_R1_001.fastq.gz, NAME_R1_001.fastq.gz)
fnFs <- sort(list.files(dir_raw, pattern="R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(dir_raw, pattern="R2_001.fastq.gz", full.names = TRUE))

# make first subdir for trimmed samples
dir_trim <- file.path(dir_cut,'trim')
if(!dir.exists(dir_trim)) dir.create(dir_trim)
# make subdir for cutadapt reports
dir_report <- file.path(dir_trim,'reports')
if(!dir.exists(dir_report)) dir.create(dir_report)

# list out file names for each sample to go into cutadapt 
fnFs_trim <- file.path(dir_trim,basename(fnFs))
fnRs_trim <- file.path(dir_trim,basename(fnRs))
# list out file names for each sample report from cutadapt
fnFRs_report <- file.path(dir_report,paste0((sapply(strsplit(basename(fnFs), "_R1_001.fastq.gz"), `[`, 1)),'_report.cutadapt.json'))

# combine F/R primers and their RC for R1 and R2 (forward and reverse reads)
R1.flags <- paste('-g',FwdPrimer,
                  '-a',RevPrimerRC)
R2.flags <- paste('-G',RevPrimer,
                  '-A',FwdPrimerRC)

# run cutadapt for each sample and write to dir_trim
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, 
                             R2.flags,
                             "--minimum-length",1, #by default, cutadapt retains all processed reads, so increasing the limit to 1 ensures downstream processing still works
                             "--discard-untrimmed",
                             "--cores",0, #automatically detects the number of available cores to enable parallel processing
                             "--json",fnFRs_report[i], #report outfile
                             "--max-n 0", #discards reads containing more than  0 N (ambiguous) bases
                             #--overlap 3, #may need to change the minimum overlap length between read and the adapter sequence. default is 3 bases
                             "-o", fnFs_trim[i], #out.1.fastq
                             "-p", fnRs_trim[i], #out.2.fastq 
                             fnFs[i], #input forward files
                             fnRs[i] #input reverse files
                             ))
}

# determine how much to trim each based on quality profile graphs
# check quality of the F of samples 3 & 4
plotQualityProfile(fnFs_trim[3:4]) 
ggsave(filename=file.path(dir_cut,"Read_F_quality_profile_trim.pdf"))
# check quality of the R of samples 3 & 4
plotQualityProfile(fnRs_trim[3:4]) 
ggsave(filename=file.path(dir_cut,"Read_R_quality_profile_trim.pdf"))


# now go run 02b_fastqc.sh to check out where to truncate reads

# Step 3: QC 
trimfnFs <- sort(list.files(dir_trim, pattern="R1_001.fastq.gz", full.names = TRUE))
trimfnRs <- sort(list.files(dir_trim, pattern="R2_001.fastq.gz", full.names = TRUE))

#make names for filtered files
sample.names <- sapply(strsplit(basename(fnFs), "_R1_001.fastq.gz"), `[`, 1)
filtFs <- file.path(dir_trim, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(dir_trim, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# filter and trim
out <- filterAndTrim(
  trimfnFs, 
  filtFs, 
  trimfnRs, 
  filtRs, 
  truncLen=c(260,260), #truncate reads after n bases
  #trimLeft=0, #number of nucleotides to remove from the start of each read. if both truncLen and trimLeft are provided, filtered reads will have length truncLen-trimLeft
  maxN=0, #DADA2 requires no Ns
  maxEE=c(2,5), #default DADA2 #max expected errors 
  truncQ=2, #default DADA2
  rm.phix=TRUE, 
  compress=TRUE, 
  multithread=TRUE, 
  verbose=TRUE) 

#preview the trimming 
head(out)
#If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads (eg. maxEE=c(2,5)), and reducing the truncLen to remove low quality tails.
#mutate(as.data.frame(out),perc=reads.out/reads.in)
#truncLen=c(275,275) less than 10% of reads are passing...
#truncLen=c(275,275),maxEE=c(2,5) less than 10% of reads are passing...
#truncLen=c(260,260),maxEE=c(2,5) >90% of reads pass

# Step 4: learn errors
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# plot errors (sanity check) 
## red line shows error rates expected under nominal def of Q score; estimated error rates are black line - ideally, estimated error rates should be a good fit for the observed error rates (points)
plotErrors(errF,nominalQ=TRUE) 
ggsave(filename=file.path(dir_trim,"ErrorRates_Plot_F.pdf"))
plotErrors(errR,nominalQ=TRUE)
ggsave(filename=file.path(dir_trim,"ErrorRates_Plot_R.pdf"))

# Step 5: denoise
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# inspect denoising
# will tell us how many sequence variants were detected from how many sequences
dadaFs[[1]] 
dadaRs[[1]]

# Step 6: merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
## most of the reads should successfully merge; if they do not we should revisit upstream parameters

# Step 7: construct sequence table 
## this returns a matrix with rows corresponding to samples, columns corresponding to sequence variants 
seqtab <- makeSequenceTable(mergers)

## check the sequence table
dim(seqtab) #rows columns
table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths

# Step 8: remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim))) # Inspect distribution of sequence lengths
## check how many chimeras/what proportion of reads they make up 
sum(seqtab.nochim)/sum(seqtab) #inspect the proportion of reads made up by chimeras


# Step 9: Track reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Step 10: save Rds files 
dir_results <- paste0(dir_cut,'/results')
if(!dir.exists(dir_results)) dir.create(dir_results)

saveRDS(out,file=paste0(dir_results,'/out_filterandtrim.rds'))
saveRDS(errF,file=paste0(dir_results,'/errF.rds'))
saveRDS(errR,file=paste0(dir_results,'/errR.rds'))
saveRDS(dadaFs,file=paste0(dir_results,'/dadaFs.rds'))
saveRDS(dadaRs,file=paste0(dir_results,'/dadaRs.rds'))
saveRDS(mergers,file=paste0(dir_results,'/mergers.rds'))
saveRDS(seqtab,file=paste0(dir_results,'/seqtab.rds'))
saveRDS(seqtab.nochim,file=paste0(dir_results,'/seqtab.nochim.rds'))
saveRDS(track,file=paste0(dir_results,'/track_reads.rds'))

## also save some files as csv
write.csv(seqtab.nochim,file=paste0(dir_results,'/seqtab.nochim.csv'))
write.csv(track,file=paste0(dir_results,'/track_reads.csv'))

## Save fasta
uniquesToFasta(getUniques(seqtab.nochim), 
               fout= paste0(dir_results,"/dada2-uniqueseqs.fasta")
               )
uniqueseqs <- readDNAStringSet(paste0(dir_results,'/dada2-uniqueseqs.fasta'))

# Step 11: return info 
Sys.Date()
sessionInfo()
save.image("/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI/BZrookery_eDNA.RData")
save.image(paste0(dir_data,"/02_Trim_QC.RData"))

# assembled with help from https://benjjneb.github.io/dada2/tutorial.html