#!/bin/bash
# Map trimmed reads to reference genome with bowtie2 
# Note: Requires that reference is indexed (see README.md) 

# Load modules
module load bowtie/2.1.0

# Global variables
TRIMMED_FOLDER="03_trimmed"
MAPPED_FOLDER="04_mapped"
REFERENCE="/home/bensuth/00_resources/sfontinalis_contigs.fasta"

# User variables
NUM_THREADS="10"

# Map reads and add RG
ls -1 $TRIMMED_FOLDER/*.fastq.gz |
	sort -u |
	while read i
	do
	  echo $i
	  name=$(basename $i)
	  label=$(echo $name | perl -pe 's/\_trimmed\.fastq\.gz//')
	  ID="@RG\tID:${label}\tSM:${label}\tPL:Illumina"
	  bowtie2 --end-to-end -k 40 --threads $NUM_THREADS --rg-id $ID -x $REFERENCE -U $i -S $i.bowtie2.sam
	  samtools view -Sb $i.sam > $i.unsorted.bam
	  samtools sort -n -o $i.sorted.bam $i.unsorted.bam
	  samtools index $i.sorted.bam
done

# clean up space
#mv ./$TRIMMED_FOLDER/*.sam ./$MAPPED_FOLDER/