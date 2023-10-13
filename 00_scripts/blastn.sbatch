#!/bin/bash
#SBATCH -J BLASTn
#SBATCH -q inferno
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1 --ntasks-per-node=2
#SBATCH --mem=15G
#SBATCH -t 24:00:00
#SBATCH -e BLASTn.err
#SBATCH -o BLASTn.out

cd $SLURM_SUBMIT_DIR

module load blast-plus/2.10.1

if [ ! -d ${odir} ]; then mkdir ${odir}; fi

n=`basename $qry | cut -d. -f1`

blastn -num_threads 2 -max_target_seqs 10 \
	-db ${ref} -query ${qry} -out ${odir}/${n}.blastn  -subject_besthit -outfmt \
	'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

scripts="/storage/scratch1/9/rconrad6/ROCkOut/00c_scripts"
python ${scripts}/besthit_filter_blast.py -i ${odir}/${n}.blastn

## Run Log
## sbatch -e BLASTn.err -o BLASTn.out --export ref=,qry=,odir= ../../00b_sbatch/blastn.sbatch
