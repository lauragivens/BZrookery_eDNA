#!/bin/bash

cd ~./blastdb/core_nt

#to create a subset of a blast database
#https://www.ncbi.nlm.nih.gov/books/NBK569848/#ckbk_aliastool.Create_a_subset_of_a_BLAS

## goal is to have a couple different subsets to BLAST off of 

# First, based on taxonomy ID 
blastdb_aliastool -db core_nt -taxidlist /Volumes/Fuji/Mangroves/NCBI_taxid_MAR.txt -dbtype nucl -out MAR_taxid 
#Created nucleotide BLAST (alias) database MAR_taxid with 865989 sequences

# Second, the mitogenomes (one whole/mostly whole mitochondrial genome per species)
blastdb_aliastool -db core_nt -seqidlist /Volumes/Fuji/Mangroves/AccNo_mitogenome_MAR.txt -dbtype nucl -out MAR_mitogenome
#Created nucleotide BLAST (alias) database MAR_mitogenome with 297 sequences

# Third, the accession numbers (representing one sequence per marker per species)
blastdb_aliastool -db core_nt -seqidlist /Volumes/Fuji/Mangroves/AccNo_partial_MAR.txt -dbtype nucl -out MAR_marker
#Created nucleotide BLAST (alias) database MAR_marker with 944 sequences


#can search one or all or these 

#31 Mar 2025