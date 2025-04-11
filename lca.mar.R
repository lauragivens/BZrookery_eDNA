#!/bin/Rscript
library(biohelper)

setwd(dir_results)
lcaPident(blast_file = 'dada2.uniques.BLAST.martaxid.tsv',
          output='lca.BLAST.martaid.csv',
          pident='no')
