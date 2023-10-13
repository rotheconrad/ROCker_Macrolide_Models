#!/bin/bash
#SBATCH -J ClustalO
#SBATCH -q inferno
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1 --ntasks-per-node=12
#SBATCH --mem=10G
#SBATCH -t 4:00:00
#SBATCH -o ClustalO.out
#SBATCH -e ClustalO.err

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Load your servers clustalo module or clustalo must be available in your PATH
module load clustal-omega/1.2.4-2r3hy4

if [ ! -s ${verified}.aln ]
  then
	n1=`basename ${verified}`
	clustalo -i $verified -o ${verified}.aln --threads 12
fi

n2=`basename ${newseqs} | cut -d_ -f2-`
clustalo -i $newseqs -o $output --profile1 ${verified}.aln --threads 12

## Run Log:
## sbatch --export verified=,newseqs=,output= ../../../00b_sbatch/seq_alignment.sbatch

