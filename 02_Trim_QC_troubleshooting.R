#Load packages
library(reticulate)
library(dada2)
library(data.table)
library(tidyverse)
library(Biostrings)

# Troubleshooting
## Cutadapt results have truncated R2 reads using original script 

# activate conda environment
use_condaenv(condaenv = 'cutadapt', required=TRUE)
# set directories
setwd('/Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI')
dir_home <- getwd()
dir_raw<-paste0(dir_home,'/fastq')
dir_cut<-paste0(dir_home,'/cutadapt')
if(!dir.exists(dir_cut)) dir.create(dir_cut)
cutadapt<-"/Users/lauragivens/miniconda3/envs/cutadapt/bin/cutadapt"

# universal variables 
# get sample names from files (assuming output follows naming convention of NAME_R1_001.fastq.gz, NAME_R1_001.fastq.gz)
fnFs <- sort(list.files(dir_raw, pattern="R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(dir_raw, pattern="R2_001.fastq.gz", full.names = TRUE))

# make first subdir for trimmed samples
dir_trim <- file.path(dir_cut,'trim-trbl')
if(!dir.exists(dir_trim)) dir.create(dir_trim)
dir_report <- file.path(dir_trim,'reports')
if(!dir.exists(dir_report)) dir.create(dir_report)

# list out file names for each sample to go into cutadapt 
fnFs_trim <- file.path(dir_trim,basename(fnFs))
fnRs_trim <- file.path(dir_trim,basename(fnRs))


# primer sequences

# sequences and citations backtracked from SERC metabarcoding digital notebook 
# citation:https://frontiersinzoology.biomedcentral.com/articles/10.1186/1742-9994-10-34
FwdPrimer=c("GGWACWGGWTGAACWGTWTAYCCYCC") #ILL-mlLCOF1
FwdPrimerRC=dada2:::rc(FwdPrimer) 
RevPrimer=c("TANACYTCNGGRTGNCCRAARAAYCA") #ILL-jgHCO2198R
RevPrimerRC=dada2:::rc(RevPrimer)

# Original Script: Remove primers  

# combine F/R primers and their RC for R1 and R2 (forward and reverse reads)
R1.flags <- paste('-g',FwdPrimer,
                  '-a',RevPrimerRC)
R2.flags <- paste('-G',FwdPrimerRC,
                  '-A',RevPrimer)

# run cutadapt and write to dir_trim
  system2(cutadapt, args = c(R1.flags, 
                             R2.flags,
                             "--discard-untrimmed",
                             "--cores",0, #automatically detects the number of available cores to enable parallel processing
                             "--json",file.path(dir_report,paste0('original.',basename(fnFs[1]),".cutadapt.json")),
                             "--max-n 0", #discards reads containing more than  0 N (ambiguous) bases
                             #--overlap 3, #may need to change the minimum overlap length between read and the adapter sequence. default is 3 bases
                             "-o", file.path(dir_trim,paste0('original.',basename(fnFs[1])))[1], #out.1.fastq
                             "-p", file.path(dir_trim,paste0('original.',basename(fnRs[1])))[1], #out.2.fastq 
                             fnFs[1], #input forward files
                             fnRs[1] #input reverse files
                             ))
# Check  
  # determine how much to trim each based on quality profile graphs
  # check quality of the F
  plotQualityProfile(file.path(dir_trim,paste0('original.',basename(fnFs[1])))[1]) 
  # check quality of the R 
  plotQualityProfile(file.path(dir_trim,paste0('original.',basename(fnRs[1])))[1]) #returns same error 
  
########################################################

# Trial 2: use Biostrings to test against all possible configurations of F and R primers 
  
  FwdPrimer=Biostrings:::DNAString("GGWACWGGWTGAACWGTWTAYCCYCC")
  FwdOrientations= c(Forward = FwdPrimer, Complement = complement(FwdPrimer), Reverse = reverse(FwdPrimer), RevComp = reverseComplement(FwdPrimer))
  FwdPrimerOrientations=sapply(FwdOrientations,toString)
  
  RevPrimer=Biostrings:::DNAString("TANACYTCNGGRTGNCCRAARAAYCA") #ILL-jgHCO2198R
  RevOrientations= c(Forward = RevPrimer, Complement = complement(RevPrimer), Reverse = reverse(RevPrimer), RevComp = reverseComplement(RevPrimer))
  RevPrimerOrientations=sapply(RevOrientations,toString)
  
  R1.flags <- paste('-g',FwdPrimerOrientations,
                    '-a',RevPrimerOrientations)
  R2.flags <- paste('-G',FwdPrimerOrientations,
                    '-A',RevPrimerOrientations)
  
  
  system2(cutadapt, args = c(R1.flags, 
                             R2.flags,
                             "--discard-untrimmed",
                             "--cores",0, #automatically detects the number of available cores to enable parallel processing
                             "--json",file.path(dir_report,paste0('biostrings.',basename(fnFs[1]),".cutadapt.json")),
                             "--max-n 0", #discards reads containing more than  0 N (ambiguous) bases
                             #--overlap 3, #may need to change the minimum overlap length between read and the adapter sequence. default is 3 bases
                             "-o", file.path(dir_trim,paste0('biostrings.',basename(fnFs[1])))[1], #out.1.fastq
                             "-p", file.path(dir_trim,paste0('biostrings.',basename(fnRs[1])))[1], #out.2.fastq 
                             fnFs[1], #input forward files
                             fnRs[1] #input reverse files
  ))
  
  # determine how much to trim each based on quality profile graphs
  # check quality of the F
  plotQualityProfile(file.path(dir_trim,paste0('biostrings.',basename(fnFs[1])))[1]) 
  # check quality of the R 
  plotQualityProfile(file.path(dir_trim,paste0('biostrings.',basename(fnRs[1])))[1]) #error
  
########################################################
  
# Trial 3: Switch F and R for R2 flag 
  R1.flags <- paste('-g',FwdPrimer,
                    '-a',RevPrimerRC)
  R2.flags <- paste('-G',RevPrimer,
                    '-A',FwdPrimerRC)
  
system2(cutadapt, args = c(R1.flags, 
                           R2.flags,
                           "--discard-untrimmed",
                           "--cores",0, #automatically detects the number of available cores to enable parallel processing
                           "--json",file.path(dir_report,paste0('revprimer.',basename(fnFs[1]),".cutadapt.json")),
                           "--max-n 0", #discards reads containing more than  0 N (ambiguous) bases
                           #--overlap 3, #may need to change the minimum overlap length between read and the adapter sequence. default is 3 bases
                           "-o", file.path(dir_trim,paste0('revprimer.',basename(fnFs[1])))[1], #out.1.fastq
                           "-p", file.path(dir_trim,paste0('revprimer.',basename(fnRs[1])))[1], #out.2.fastq 
                           fnFs[1], #input forward files
                           fnRs[1] #input reverse files
))


# determine how much to trim each based on quality profile graphs
# check quality of the F
plotQualityProfile(file.path(dir_trim,paste0('revprimer.',basename(fnFs[1])))[1]) 
# check quality of the R 
plotQualityProfile(file.path(dir_trim,paste0('revprimer.',basename(fnRs[1])))[1]) #still error, although the number of kept reads is a lot better


########################################################

# Trial 4: uneven reads 
## uses switched F/R for R2 flag and loosen pairing filter
system2(cutadapt, args = c(R1.flags, 
                           R2.flags,
                           "--pair-filter=both", #filtering criteria must apply to both reads 
                           "--discard-untrimmed", #because pair-filter=both, the pair is discarded if both reads fail the filter
                           "--cores",0, #automatically detects the number of available cores to enable parallel processing
                           "--json",file.path(dir_report,paste0('pairboth.',basename(fnFs[1]),".cutadapt.json")),
                           "--max-n 0", #discards reads containing more than  0 N (ambiguous) bases
                           #--overlap 3, #may need to change the minimum overlap length between read and the adapter sequence. default is 3 bases
                           "-o", file.path(dir_trim,paste0('pairboth.',basename(fnFs[1])))[1], #out.1.fastq
                           "-p", file.path(dir_trim,paste0('pairboth.',basename(fnRs[1])))[1], #out.2.fastq 
                           fnFs[1], #input forward files
                           fnRs[1] #input reverse files
))
## results in a higher percentage of total written reads but still uneven between read one and read 2
plotQualityProfile(file.path(dir_trim,paste0('pairboth.',basename(fnFs[1])))[1]) 
# check quality of the R 
plotQualityProfile(file.path(dir_trim,paste0('pairboth.',basename(fnRs[1])))[1]) #still error, although the number of kept reads is a lot better

########################################################

# Trial 5: minimum length
system2(cutadapt, args = c(R1.flags, 
                           R2.flags,
                           "--minimum-length",1, #by default, cutadapt retains all processed reads, so increasing the limit to 1 ensures downstream processing still works
                           "--discard-untrimmed", #because pair-filter=both, the pair is discarded if both reads fail the filter
                           "--cores",0, #automatically detects the number of available cores to enable parallel processing
                           "--json",file.path(dir_report,paste0('minlength.',basename(fnFs[1]),".cutadapt.json")),
                           "--max-n 0", #discards reads containing more than  0 N (ambiguous) bases
                           #--overlap 3, #may need to change the minimum overlap length between read and the adapter sequence. default is 3 bases
                           "-o", file.path(dir_trim,paste0('minlength.',basename(fnFs[1])))[1], #out.1.fastq
                           "-p", file.path(dir_trim,paste0('minlength.',basename(fnRs[1])))[1], #out.2.fastq 
                           fnFs[1], #input forward files
                           fnRs[1] #input reverse files
))
## results in a higher percentage of total written reads but still uneven between read one and read 2
plotQualityProfile(file.path(dir_trim,paste0('minlength.',basename(fnFs[1])))[1]) 
# check quality of the R 
plotQualityProfile(file.path(dir_trim,paste0('minlength.',basename(fnRs[1])))[1]) #finally runs


########################################################
# FASTQC 
## outside of R environment
# conda activate fastqc
# cd /Volumes/Fuji/Mangroves/2025*
# find cutadapt/trim-trbl -type f -name "*.fastq" -o -name "*.fastq.gz" | xargs fastqc --threads 4 --outdir "cutadapt/trim-trbl/reports"
# cd cutadapt/trim-trbl/reports
# multiqc .
# multiqc *R2* --filename "cutadapt.trblshooting.R2"

## from the looks of this, pairboth has the highest number of unique reads, which is unsurprising considering the stringency to throw out pairs is much lower
## also had the lowest number of overrepresented sequences 
## minlength doesn't look like it has as much adapter contamination as the other filters; original and revprimer have high levels of adapter contamination
## minlength failed per base sequence quality and GC content, but was middling for seq length distribution
## pairboth and biostrings had better per base sequence quality

########################################################
# Trial 6: biostrings + revprimer + minimum length
## combine the best of the attempts
## reversing F/R tags in R2 seems like the way to go to retain seqs
## using biostrings orientations caught more 
R1.flags <- paste('-g',FwdPrimerOrientations,
                  '-a',RevPrimerOrientations)
R2.flags <- paste('-G',RevPrimerOrientations,
                  '-A',FwdPrimerOrientations)

system2(cutadapt, args = c(R1.flags, 
                           R2.flags,
                           "--minimum-length",1, #by default, cutadapt retains all processed reads, so increasing the limit to 1 ensures downstream processing still works
                           "--discard-untrimmed", #because pair-filter=both, the pair is discarded if both reads fail the filter
                           "--cores",0, #automatically detects the number of available cores to enable parallel processing
                           "--json",file.path(dir_report,paste0('length.biostrings.rev.',basename(fnFs[1]),".cutadapt.json")),
                           "--max-n 0", #discards reads containing more than  0 N (ambiguous) bases
                           #--overlap 3, #may need to change the minimum overlap length between read and the adapter sequence. default is 3 bases
                           "-o", file.path(dir_trim,paste0('length.biostrings.rev.',basename(fnFs[1])))[1], #out.1.fastq
                           "-p", file.path(dir_trim,paste0('length.biostrings.rev.',basename(fnRs[1])))[1], #out.2.fastq 
                           fnFs[1], #input forward files
                           fnRs[1] #input reverse files
))

plotQualityProfile(file.path(dir_trim,paste0('length.biostrings.rev.',basename(fnFs[1])))[1]) 
# check quality of the R 
plotQualityProfile(file.path(dir_trim,paste0('length.biostrings.rev.',basename(fnRs[1])))[1])

########################################################
# FASTQC 
## outside of R environment
# find cutadapt/trim-trbl -type f -name "length.biostrings.rev.*.fastq" -o -name "length.biostrings.rev.*.fastq.gz" | xargs fastqc --threads 4 --outdir "cutadapt/trim-trbl/reports"
# multiqc *R2* --filename "cutadapt.trblshooting.R2"

## doing this does improve per base sequence quality but not GC content; has slightly less dups than minlength alone,

########################################################
# Trial 7: go back to original but retain minimum length req
R1.flags <- paste('-g',FwdPrimer,
                  '-a',RevPrimerRC)
R2.flags <- paste('-G',FwdPrimerRC,
                  '-A',RevPrimer)

# run cutadapt and write to dir_trim
system2(cutadapt, args = c(R1.flags, 
                           R2.flags,
                           "--discard-untrimmed",
                           "-m",1,
                           "--cores",0, #automatically detects the number of available cores to enable parallel processing
                           "--json",file.path(dir_report,paste0('original.m..',basename(fnFs[1]),".cutadapt.json")),
                           "--max-n 0", #discards reads containing more than  0 N (ambiguous) bases
                           #--overlap 3, #may need to change the minimum overlap length between read and the adapter sequence. default is 3 bases
                           "-o", file.path(dir_trim,paste0('original.m..',basename(fnFs[1])))[1], #out.1.fastq
                           "-p", file.path(dir_trim,paste0('original.m..',basename(fnRs[1])))[1], #out.2.fastq 
                           fnFs[1], #input forward files
                           fnRs[1] #input reverse files
))
# Check  
# determine how much to trim each based on quality profile graphs
# check quality of the F
plotQualityProfile(file.path(dir_trim,paste0('original.m..',basename(fnFs[1])))[1]) 
# check quality of the R 
plotQualityProfile(file.path(dir_trim,paste0('original.m..',basename(fnRs[1])))[1]) 

## ok just wanted to double check - but yes, absolutely, this still looks like shit, retains very few reads, bad QC




# Summary
## multiqc in conda env to check all the versions that DIDNT result in ultra-truncated reads
# multiqc length*R2* *length*R2* pairboth*R2* rev*R2* --filename "multiqc.cutadapt.trblshooting.R2.length.pair.rev"
## looks at reverse primer, pair both + reverse primer, min length + reverse primer, and min length + reverse primer + biostring (not in order of above)
## and in general, it looks like the min+rev+biost is best? 
# lowest % of duplicated seqs other than pair+rev, which again we're assuming is because of the looser filtering 
# end of sequence QC drops off the slowest
# however, does have higher percentage of adapter contam, probably due to the biostrings addition 

# Final trial
########################################################
# Trial 8:  revprimer + minimum length
R1.flags <- paste('-g',FwdPrimer,
                  '-a',RevPrimerRC)
R2.flags <- paste('-G',RevPrimer,
                  '-A',FwdPrimerRC)

system2(cutadapt, args = c(R1.flags, 
                           R2.flags,
                           "--minimum-length",1, #by default, cutadapt retains all processed reads, so increasing the limit to 1 ensures downstream processing still works
                           "--discard-untrimmed", #because pair-filter=both, the pair is discarded if both reads fail the filter
                           "--cores",0, #automatically detects the number of available cores to enable parallel processing
                           "--json",file.path(dir_report,paste0('length.rev.',basename(fnFs[1]),".cutadapt.json")),
                           "--max-n 0", #discards reads containing more than  0 N (ambiguous) bases
                           #--overlap 3, #may need to change the minimum overlap length between read and the adapter sequence. default is 3 bases
                           "-o", file.path(dir_trim,paste0('length.rev.',basename(fnFs[1])))[1], #out.1.fastq
                           "-p", file.path(dir_trim,paste0('length.rev.',basename(fnRs[1])))[1], #out.2.fastq 
                           fnFs[1], #input forward files
                           fnRs[1] #input reverse files
))

plotQualityProfile(file.path(dir_trim,paste0('length.rev.',basename(fnFs[1])))[1]) 
# check quality of the R 
plotQualityProfile(file.path(dir_trim,paste0('length.rev.',basename(fnRs[1])))[1])

########################################################
# FASTQC 
## outside of R environment
# find cutadapt/trim-trbl -type f -name "length.rev.*.fastq" -o -name "length.rev.*.fastq.gz" | xargs fastqc --threads 4 --outdir "cutadapt/trim-trbl/reports"
# multiqc length*R2* *length*R2*  rev*R2* --filename "multiqc.cutadapt.trblshooting.R2.length.rev"

# Summary
## mainly interested in comparing the rev+minlength+biostrings and rev+minlength, but leaving the separate ones in for context
## we see that length+biostrings+rev has lowest % duplicated reads (barely), all have ~same number of reads
## last ~20bp or so on all four versions is pretty shit
## biostrings version still has the highest  adapter contamination, but all pass 

# So I don't really see a pressing reason to keep the Biostrings version at this exact moment. Will proceed with the rev+minlength

Sys.Date()
sessionInfo()
