#!/bin/bash
#SBATCH --job-name=filtering
#SBATCH -n 24
#SBATCH --partition=compute
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --output=filter.out
#SBATCH --error=filter.err

ml samtools

for bam in *Aligned.out.bam
do
    # Create output filename
    base=$(basename "$bam" .bam)
    out="${base}_mapq20.bam"

    echo "Filtering $bam -> $out (keeping MAPQ >= 20)"

    # Filter and write to new BAM
    samtools view -@ 8 -b -q 20 "$bam" -o "$out"

    # Index the filtered BAM
    samtools index "$out"
done
