#!/bin/bash
#SBATCH -J IQTree
#SBATCH -q inferno
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1 --ntasks-per-node=20
#SBATCH --mem=50G
#SBATCH -t 5-00:00:00
#SBATCH -o IQTree.out
#SBATCH -e IQTree.err

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Load your servers IQTREE module or make IQTREE available in your PATH
# IQTREE works well with a conda install
# load conda env for rockout
source ~/.bashrc
source activate /storage/home/hcoda1/9/rconrad6/p-ktk3-0/apps/ROCkOut

iqtree -s ${input} --prefix ${outpre} -m MFP -alrt 1000 -B 1000 -T 20

## Run Log ##
## sbatch --export input=fasta.aln,outpre= /Path/to/GitHub/repo/01b_Sbatch/03c_IQTree.sbatch

