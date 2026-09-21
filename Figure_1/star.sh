#!/bin/bash
#SBATCH --job-name=STAR
#SBATCH --cpus-per-task=12
#SBATCH --partition=gen-mk-compute-1
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --array=1-66%5
#SBATCH --output=logs/star_%A_%a.out
#SBATCH --error=logs/star_%A_%a.err

# Load software
module load star/2.7.10a
ml samtools

name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" complete_samples.txt) #complete_samples.txt is list of abbreviated file names

STAR --runThreadN 12 --runMode alignReads \
--twopassMode Basic \
--outSAMtype BAM Unsorted \
--readFilesCommand gunzip -c \
--genomeDir /data2/lackey_lab/DownloadedSequenceData/randazza/new_chick/star/ \
--outFileNamePrefix ${name}_ \
--readFilesIn ${name}_clean_R1.fastq.gz ${name}_clean_R2.fastq.gz \
--outFilterType BySJout \
--outSAMattributes NH HI AS NM MD jM jI \
--outFilterMultimapNmax 100 \
--winAnchorMultimapNmax 200 \
--outFilterMismatchNoverReadLmax 0.04 \
--outReadsUnmapped Fastx \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--quantMode TranscriptomeSAM \
--alignSJDBoverhangMin 1
