#!/bin/bash
  
#SBATCH --job-name=rookeryBLASTalias
#SBATCH --account=schultzlab -p schultzlab
#SBATCH --output=/work/lag66/rookery/alias_ntBLAST.out
#SBATCH --error=/work/lag66/rookery/alias_ntBLAST.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=150G
#SBATCH --mail-user=lag66@duke.edu

# --------------------------------------- load modules ---------------------------------------
module load NCBI-BLAST/2.12.0-rhel8
echo loaded BLAST module 

FILES="/work/lag66/rookery"

# --------------------------------------- set path ---------------------------------------
source /hpc/group/schultzlab/lag66/miniconda3/etc/profile.d/conda.sh

# --------------------------------------- assign directories ---------------------------------------
cd /hpc/group/schultzlab/lag66
export PATH="$PATH:/hpc/group/schultzlab/lag66/core_nt"
BLASTDB="/hpc/group/schultzlab/lag66/core_nt"
export BLASTDB
#TAXDB="/hpc/group/schultzlab/lag66"
 
echo we will be BLASTing against a local BLAST nt database


# --------------------------------------- start BLAST ---------------------------------------
echo starting BLAST against MAR alias db  
echo 
blastn -db core_nt/MAR_taxid -query $FILES/dada2-uniqueseqs.fasta \
                -perc_identity 97 \
                        -word_size 30 \
                                -outfmt 6 \
                                        -num_threads 4 \
                                                -out $FILES/dada2.uniques.BLAST.martaxid.tsv


echo $(date)
