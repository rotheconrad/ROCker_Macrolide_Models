#!/bin/bash
#SBATCH -J SixFrames
#SBATCH -q inferno
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH -t 12:00:00
#SBATCH -e SixFrames.err
#SBATCH -o SixFrames.out

cd $SLURM_SUBMIT_DIR

script="/storage/scratch1/9/rconrad6/ROCkOut/00c_scripts/nuc_fasta_6aa_frame_fasta.py"

python ${script} -i ${infile} -o ${outfile}

## Run Log
## sbatch --export infile=,outfile= ../../00b_sbatch/TranslateSixFrames.sbatch
## sbatch --export infile=mock_metagenome_labeled.fasta,outfile=mock_metagenome_labeled_6frames.faa ../../00b_sbatch/FragGeneScanRs.sbatch
