#!/bin/bash
#SBATCH --job-name=STAR
#SBATCH -n 24
#SBATCH --partition=compute
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --output=/data2/lackey_lab/DownloadedSequenceData/randazza/chick_seq/merged_nova/star.out
#SBATCH --error=/data2/lackey_lab/DownloadedSequenceData/randazza/chick_seq/merged_nova/star.err

# Load software
cd /data2/lackey_lab/DownloadedSequenceData/randazza/chick_seq/merged_nova/lower_map_align/single_run 
module load star/2.7.10a
ml samtools

for i in *_R1_001.fastq.gz; do name=$(basename ${i} _R1_001.fastq.gz);
STAR --runThreadN 24 --runMode alignReads \
--twopassMode Basic \
--outSAMtype BAM Unsorted \
--readFilesCommand gunzip -c \
--genomeDir /data2/lackey_lab/DownloadedSequenceData/randazza/new_chick/star/ \
--outFileNamePrefix ${name}_ \
--readFilesIn ${name}_R1_001.fastq.gz ${name}_R2_001.fastq.gz \
--outFilterType BySJout \
--outSAMattributes NH HI AS NM MD jM jI \
--outFilterMultimapNmax 20 \
--outFilterMismatchNoverReadLmax 0.04 \
--outReadsUnmapped Fastx \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--quantMode TranscriptomeSAM \
--outFilterScoreMinOverLread 0.3 \
--outFilterMatchNminOverLread 0.3 \
--alignSJDBoverhangMin 1 ;
done
