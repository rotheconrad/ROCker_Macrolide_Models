#!/bin/bash
#SBATCH -J BLASTx
#SBATCH -q inferno
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1 --ntasks-per-node=2
#SBATCH --mem=15G
#SBATCH -t 24:00:00
#SBATCH -e BLASTx.err
#SBATCH -o BLASTx.out

cd $SLURM_SUBMIT_DIR

module load blast-plus/2.10.1

if [ ! -d ${odir} ]; then mkdir ${odir}; fi

n=`basename $qry | cut -d. -f1`

blastx -num_threads 2 -max_target_seqs 10 \
	-db ${ref} -query ${qry} -out ${odir}/${n}.blastx  -subject_besthit -outfmt \
	'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

scripts="/storage/scratch1/9/rconrad6/ROCkOut/00c_scripts"
python ${scripts}/besthit_filter_blast.py -i ${odir}/${n}.blastx

## Run Log
## sbatch -e BLASTx.err -o BLASTx.out --export ref=,qry=,odir= ../../00b_sbatch/blastx.sbatch
