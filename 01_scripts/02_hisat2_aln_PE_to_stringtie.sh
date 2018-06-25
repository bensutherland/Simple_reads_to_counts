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
REFERENCE="/hdd/genomes/GCF_002872995.1_Otsh_v1.0_genomic.fna"

# User variables
NUM_THREADS="8"

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
	  hisat2 -x $REFERENCE -k 40 --threads $NUM_THREADS --rg-id $ID -1 $i"_1.paired.fastq.gz" -2 $i"_2.paired.fastq.gz" -S $i.hisat2.sam
	  samtools view -Sb $i.hisat2.sam > $i.hisat2.unsorted.bam
	  samtools sort -o $i.hisat2.sorted.bam $i.hisat2.unsorted.bam
	  #samtools index $i.bowtie2.sorted.bam # note, if samtools sort -n is used, cannot index
      rm $i.hisat2.sam
done

# Move output bam files into the mapped folder
mv $TRIMMED_FOLDER/*.bam $MAPPED_FOLDER
