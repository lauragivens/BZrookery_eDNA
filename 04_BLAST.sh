#!/bin/bash 

# updated script to use remote tag for BLAST

############# Download and install command line BLAST+ ############# 
# Windows installer: https://www.ncbi.nlm.nih.gov/books/NBK52637/
# LINUX/UNIX installer: https://www.ncbi.nlm.nih.gov/books/NBK52640/
# Mac installer: https://www.ncbi.nlm.nih.gov/books/NBK569861/

# dmg files are installers
# tar files require compilation 

# after installation, you will more than likely have to 
# configure your machine by adding BLAST to your path, 
# aka where your machine looks for the package.  
# https://www.ncbi.nlm.nih.gov/books/NBK52640/

# switch the BLAST name to whichever version you downloaded  
export PATH=$PATH:$HOME/ncbi-blast-2.17.0+/bin

# then you can check that it was added with 
echo $PATH

# next, add the location where you will download any BLAST databases 
export BLASTDB=$HOME/blastdb

# Notice above that the BLAST package is being pointed to $HOME/blast/bin
# And your databases are going into $HOME/blastdb  
# these are not the same things 
# the package is what your machine will use to perform the BLAST search
# the databases are what you compare your sequences against

# check that setting your path works  
blastn -help 
# if your machine does not recognize the command, setting the path did not work

############# Download databases ############# 
# this is required to run BLAST locally
# at minimum you should have one of the BLAST databases downloaded
# and the taxonoy database 
# https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz

# this may take a lot of memory 

### WARNING ###
# DO NOT DO THIS FOR THIS TUTORIAL # 
# SKIP TO REMOTE SECTION # 

# Download a pre-formatted database 
#update_blastdb.pl --decompress core_nt 
# Download the taxdb archive
#update_blastdb.pl taxdb
# Install it in the BLASTDB directory
#gunzip -cd taxdb.tar.gz | (cd $BLASTDB; tar xvf - )

############# Remote option ############# 
# add -remote tag to the end of your query 
# more info here: 
# https://www.ncbi.nlm.nih.gov/books/NBK569854/#ckbk_DELTABLASTsrch.BLAST_remote_service

############# Taxonomic assignment #############  
# quick start: https://www.ncbi.nlm.nih.gov/books/NBK569856/
# basic structure is
# blastn -db nt -query in.fasta -out results.out
FILES="/Users/lauragivens/Downloads/Demo/cutadapt/results"

blastn -db core_nt -query $FILES/dada2-uniqueseqs.fasta \
-perc_identity 97 \
-qcov_hsp_perc 80 \
-word_size 30 \
-outfmt 6 \
-out $FILES/dada2.uniques.BLAST.default.tsv \
-remote

# this will then run and save the output to the file you defined as out
