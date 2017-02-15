#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH -J samtools 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bensutherland7@gmail.com
#SBATCH --time=44:00:00
#SBATCH --mem=10000

# Filter mapped reads and make sorted, indexed bam file

# Load modules
module load samtools/1.3

# Global variables
MAPPED_FOLDER="04_mapped"

# User variables
NUM_THREADS="2"


# map reads and add RG
ls -1 $MAPPED_FOLDER/*.sam |
	sort -u |
	while read i
	do
	  echo $i
	  samtools view -Sb $i > $i.unsorted.bam
	  samtools sort -o $i.sorted.bam $i.unsorted.bam
	  samtools index $i.sorted.bam
done
