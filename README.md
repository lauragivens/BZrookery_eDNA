# BZrookery_eDNA
analysis of rookery/non-rookery mangrove eDNA samples 

This is a workflow outlining the bioinformatic processing steps taken after sequencing results were returned. 

You will need the following: 
## Conda environment  
https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html  
 
### fastqc  
Download instructions: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt

### multiqc 
Github: https://github.com/MultiQC/MultiQC  

### cutadapt  
Github: https://github.com/marcelm/cutadapt/tree/main

## R packages  
reticulate
dada2
data.table
tidyverse
Biostrings
lulu
taxonomizr  
phyloseq  
vegan  
decontam  
ggplot2  
cowplot  

## Command Line BLAST  
https://www.ncbi.nlm.nih.gov/books/NBK569861/  


# 01_fastqc.sh  
In this script, we use a conda environment to activate *FASTQC* and check every .fastq and zipped .fastq.gz file from the sequencing center.  
Then, we run *multiqc* on the reports to aggregate and summarize the results.  

You can find the aggregated report in the FASTQC folder.  

# 02a_Trim_QC.R  
This script calls the cutadapt conda environment from within R, which requires you to tell R to use a conda environment and then point to the location of that conda environment on your machine.  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use_condaenv(condaenv = 'cutadapt', required=TRUE)
cutadapt<-"/Users/lauragivens/miniconda3/envs/cutadapt/bin/cutadapt"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You can then use the environment.   
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
system2(cutadapt, args = c(FwdPrimer,RevPrimer,
FwdReads,RevReads,
-o out.f.fastq,
-p out.r.fastq
))
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
It takes the foward and reverse primer sequences and their reverse complements and uses cutadapt to trim them from each read.  
It then will plot the quality profile of the reads, which is then used to determine how much to truncate the reads. You want to truncate enough that you are not including too many low-quality bases, but need to retain enough overlap that you can merge the forward and reverse reads later.   

# 02b_fastqc.R
We once again go back to fastqc and check the quality of the forward and reverse reads. Using the information from fastqc and the quality profile plot from cutadapt, we can then define how many bases to trim from each end of the reads.  

# 02a_Trim_QC.R  
Reads are then filtered and trimmed based on user defined standards.
Reads are denoised and merged, then we make a sequence table and remove chimeras.  

# 03_Cluster.R   
 This script uses LULU (https://github.com/tobiasgf/lulu) to cluster sequences into ASVs 
This used BLASTn to produce a match list, which blasts the sequences against each other. 

# 04_runBLAST**.sh   
This was run on a HPC  
BLAST aliases were created using the Create-subdb.sh script in the SetUp_Scripts folder and the id lists in the SetUp_Files folder.  
This script can be run on a local machine with some edits, but may run into memory problems.  

# 05_AccNo_toTax.R   
The output of the previous step does not include the taxonomic ranks. This script uses the taxonomizr R package (https://cran.r-project.org/web/packages/taxonomizr/readme/README.html) to convert accession numbers to taxonomic ranks  

# 06_toPhyloseq.R  
First, this script will truncate taxonomic assignments based on percent match value from BLAST.  
It then adds trophic information and other background data to the taxonomy table based on taxid from NCBI.  
This script is used to upload everything into phyloseq for visualization.  

# 07_Blanks.Rmd  
After data are uploaded into phyloseq format, we then use the decontam R package (https://github.com/benjjneb/decontam) to identify contaminants based on the distribution of sequence frequency and prevalence in the various control blanks.  
After identifying and removing contaminants, we then identify each ASV in each blank and its abundance within the blank. These ASVs are then removed from the corresponding samples.  
For example, for an ASV identified in a field blank for a specific sampling day, the abundance of that ASV within that blank is removed from each sample collected on the same sampling day.  
Thus, the abundance of each sequence occurring in a field or lab blank is removed from all samples collected on that sampling day. The abundance of each sequence occurring in an extraction blank was removed from all samples extracted on the same day. The abundance of each sequence within any PCR blank was removed from all samples.      

