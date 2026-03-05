#!/bin/bash
#SBATCH --job-name=partition
#SBATCH -n 16
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --error=partition.err
# Load software
ml rnastructure
export DATAPATH=/data/software/RNAstructure_v6.5/data_tables/

cd Aexon/
for seq in *.fasta; do
	N=$(basename ${seq} .fasta);
	partition $seq ../microA3_bp/$N.pfs -md 66 
done
