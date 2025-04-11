library(biohelper)

#wrapper function to perform a blastn search then uses a lca approach and optionally a minimum percent identity to assign taxonomy
blastn_taxo_assignment(blastapp_path = '/Users/lauragivens/ncbi-blast-2.16.0+/bin', #path to blastn program
                       db='core_nt/MAR_taxid',#reference database
                       queries=paste0(dir_results,'/dada2-uniqueseqs.fasta'),
                       method='blastn', #default='both'
                       output_path=paste0(dir_results,'/blast_results'),
                       nthreads=10,
                       pident='no' #default
                       )
