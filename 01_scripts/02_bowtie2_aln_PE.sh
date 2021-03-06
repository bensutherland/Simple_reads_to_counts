#!/bin/bash
# Map paired-end trimmed reads to reference genome with bowtie2 
# Note: Requires that reference is indexed (see README.md) 
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Global variables
TRIMMED_FOLDER="03_trimmed"
MAPPED_FOLDER="04_mapped"
REFERENCE="/home/ben/Documents/genomes/GCF_002872995.1_Otsh_v1.0_genomic.fna"

# User variables
NUM_THREADS="6"

# Map reads and add RG
ls -1 $TRIMMED_FOLDER/*.paired.fastq.gz |
	perl -pe 's/\_[12]\.paired\.fastq\.gz//' |
	sort -u |
	while read i
	do
	  echo $i
	  name=$(basename $i)
	  label=$(echo $name)
	  ID="@RG\tID:${label}\tSM:${label}\tPL:Illumina"
	  bowtie2 --end-to-end -k 40 --threads $NUM_THREADS --rg-id $ID -x $REFERENCE -1 $i"_1.paired.fastq.gz" -2 $i"_2.paired.fastq.gz" -S $i.bowtie2.sam
	  samtools view -Sb $i.bowtie2.sam > $i.bowtie2.unsorted.bam
	  samtools sort -n -o $i.bowtie2.sorted.bam $i.bowtie2.unsorted.bam
	  #samtools index $i.bowtie2.sorted.bam # note, if samtools sort -n is used, cannot index
      rm $i.bowtie2.sam
done

# Move output bam files into the mapped folder
mv $TRIMMED_FOLDER/*.bam $MAPPED_FOLDER
