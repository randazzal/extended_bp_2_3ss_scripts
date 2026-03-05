#!/bin/bash
#SBATCH --job-name=partition
#SBATCH -n 16
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --error=pairprob.err
# Load software
ml rnastructure
export DATAPATH=/data/software/RNAstructure_v6.5/data_tables
cd microA3_bp/

for file in *.pfs; do
	N=$(basename ${file} .pfs);
	ProbabilityPlot $file $N.dp -t
done
