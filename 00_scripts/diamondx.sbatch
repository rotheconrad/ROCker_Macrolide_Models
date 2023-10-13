#!/bin/bash
#SBATCH -J DIAMONDx
#SBATCH -q inferno
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1 --ntasks-per-node=2
#SBATCH --mem=15G
#SBATCH -t 24:00:00
#SBATCH -e DIAMONDx.err
#SBATCH -o DIAMONDx.out

cd $SLURM_SUBMIT_DIR

# load conda env for rockout
source ~/.bashrc
source activate /storage/home/hcoda1/9/rconrad6/p-ktk3-0/apps/ROCkOut

echo ' '
echo Conda env: $CONDA_PREFIX

if [ ! -d ${odir} ]; then mkdir ${odir}; fi

n=`basename $qry | cut -d. -f1`

diamond blastx --very-sensitive --unal 0 --threads 2 \
	--db ${ref} --query ${qry} --out ${odir}/${n}.dmndx \
	--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

if [ -s ${odir}/${n}.dmndx ]; then
	scripts="/storage/scratch1/9/rconrad6/ROCkOut/00c_scripts"
	python ${scripts}/besthit_filter_blast.py -i ${odir}/${n}.dmndx
  fi

## Run Log
## sbatch -e DIAMONDx.err -o DIAMONDx.out --export ref=,qry=,odir= ../../00b_sbatch/diamondx.sbatch
