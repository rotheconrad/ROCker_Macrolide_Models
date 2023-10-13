#!/bin/bash
#SBATCH -J alignTests
#SBATCH -q inferno
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1 --ntasks-per-node=2
#SBATCH --mem=15G
#SBATCH -t 24:00:00
#SBATCH -e alignment_label_test.err
#SBATCH -o alignment_label_test.out

cd $SLURM_SUBMIT_DIR

script="/storage/scratch1/9/rconrad6/ROCkOut/00c_scripts/nuc_fasta_6aa_frame_fasta.py"

module load blast-plus/2.10.1

if [ ! -d ${odirn} ]; then mkdir ${odirn}; fi
if [ ! -d ${odira} ]; then mkdir ${odira}; fi

for f in 01a_self_score_fasta/*.fasta;
	do
		n=`basename $f | cut -d. -f1`;
		blastn -num_threads 2 -max_target_seqs 10 \
		-db ${refn} -query ${f} -out ${odirn}/${n}.blastn  -subject_besthit -outfmt \
		'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

		blastx -num_threads 2 -max_target_seqs 10 \
		-db ${refa} -query ${f} -out ${odira}/${n}.blastx  -subject_besthit -outfmt \
		'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

	done


#python ${script} -i ${infile} -o ${outfile}

## Run Log
## sbatch -e self_score.err -o self_score.out --export refn=,refa=,fdir=,odirn=,odira= ../../00b_sbatch/alignment_label_tests.sbatch
