#!/bin/bash
#SBATCH -J TrimAl
#SBATCH -q inferno
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem=1G
#SBATCH -t 12:00:00
#SBATCH -o seq_trim.out
#SBATCH -e seq_trim.err

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Load your servers trimal module, make trimal available in your path, or provide the path to trimal as a variable "trimal"
trimal=/storage/coda1/p-ktk3/0/rconrad6/apps/trimal/source/./trimal

${trimal} -in ${input} -out ${output} -resoverlap 0.75 -seqoverlap 80 -automated1

## Param key:
## -resoverlap: Minimum overlap of a positions with other positions in the column to be considered a "good position". Range: [0 - 1]. (see User Guide).
## -seqoverlap: Minimum percentage of "good positions" that a sequence must have in order to be conserved. Range: [0 - 100](see User Guide).

## Run Log ##
## sbatch --export input=fasta.aln,output=fasta.aln /Path/to/GitHub/repo/01b_Sbatch/02e_seq_trim.sbatch
## sbatch --export input=final_tree_seqs.aln,output=final_tree_seqs_trimmed.aln ../../../00b_sbatch/seq_trim.sbatch
