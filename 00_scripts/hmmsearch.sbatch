#!/bin/bash
#SBATCH -J hmmsearch
#SBATCH -q inferno
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1 --ntasks-per-node=20
#SBATCH --mem=50G
#SBATCH -t 24:00:00
#SBATCH -e hmmsearch.err
#SBATCH -o hmmsearch.out

cd $SLURM_SUBMIT_DIR

# load conda env for rockout
source ~/.bashrc
source activate /storage/home/hcoda1/9/rconrad6/p-ktk3-0/apps/ROCkOut

echo ' '
echo Conda env: $CONDA_PREFIX

hmmsearch --cpu 20 ${model} ${infile} > ${outfile}

## Run Log
## sbatch --export model=,infile=,outfile= ../../00b_sbatch/hmmsearch.sbatch
## sbatch --export model=mcrAll_model.hmm,infile=mock_metagenome_labeled_6frames.faa,outfile=mcrALL_hmmSearch.tsv ../../00b_sbatch/hmmsearch.sbatch
