#!/bin/bash 

#https://github.com/DominikBuchner/BOLDigger3

conda create -n py3129 python=3.12.9
conda activate py3129
pip install boldigger3
playwright install

cd /Volumes/Fuji/Mangroves/2025_0319_Givens_Canty_Rookery_COI/cutadapt/results

#basic structure
#boldigger3 identify path_to_pasta --db database_nr --mode operating_mode
#db refers to which bold v5 db to use 
#operating_mode refers to whether to use rapid search, genus-species search, or exhaustive

boldigger3 identify dada2-uniqueseqs.fasta --db 2 --mode 2