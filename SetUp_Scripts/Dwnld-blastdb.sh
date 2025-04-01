#!/bin/bash

#steps to download ncbi blast package mar 27 2025
#https://www.ncbi.nlm.nih.gov/books/NBK52640/
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-macosx.tar.gz

#To install, simply extract the downloaded package after placing it under a desired directory. This can be accomplished by a single tar command, or a combination of gunzip and tar commands
tar zxvpf ncbi-blast-2.16.0+-x64-macosx.tar.gz

#Successful execution of the above commands installs the package and generates a new ncbi-blast-2.10.1+ directory under the working directory selected. This new directory contains the bin and doc subdirectories, as well as a set of informational files.

#Using the BLAST+ package installed above without configuration will be cumbersome – it requires the installation path to be prefixed to the program and database calls since the system does not know where to look for the installed program and the specified database. To streamline BLAST searches, two environment variables, PATH and BLASTDB, need to be modified and created, respectively, to point to the corresponding directories.
export PATH=$PATH:$HOME/ncbi-blast-2.16.0+/bin

#create a directory to store BLAST databases
mkdir $HOME/blastdb
export BLASTDB=$HOME/blastdb

#Once they are set, the system knows where to call BLAST programs, and the invoked program will know where to look for the database files. 
#
###### BLAST DB download
#https://www.ncbi.nlm.nih.gov/books/NBK569850/
#
cd $HOME/blastdb

#now to download the nt database
#https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html
#Pre-formatted databases must be downloaded using the update_blastdb.pl script or via FTP in binary mode.

#check available databases for download
update_blastdb.pl --showall [*] #this will connect to ncbi 
# possible options include core_nt, nt_euk, or mito? 

#based on blasting one mito accession no for core and euk, I think nt_euk is the smallest. mito isn't an option for web BLAST
#update_blastdb.pl --decompress nt_euk [*]
#ok that's way too fucking big; 181 files at least 2.5 G each
#update_blastdb.pl -decompress core_nt [*]
update_blastdb.pl core_nt [*]

#does not decompress and throws error
#use this from when I downloaded nt db to DCC
 ls *.gz |xargs -n1 tar -xzvf
 mkdir core_nt_gz
  mv *.gz core_nt_gz #originally deleted these files, but it seems like update_blastdb.pl needs the original tar.gz files to check timestamps to determine if a refresh is needed
  
  #check database downloaded and unzipped correctly
  blastdbcheck -db core_nt 
  
  # move to subdirectory
  mkdir core_nt
  mv *.* core_nt
  
  #i don't THINK i need to download the tax db separately because of the following from the cookbook: 
  #If you do add the taxids to your database, make sure you have the BLAST taxonomy data files (taxdb.bt[di]) which are available from https://ftp.ncbi.nlm.nih.gov/blast/db/ but also packaged with most BLAST databases distributed by the NCBI.
  
  ## Downloaded 31 Mar 2025