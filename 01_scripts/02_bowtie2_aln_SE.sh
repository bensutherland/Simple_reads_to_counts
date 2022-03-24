#!/bin/bash
# Map trimmed reads to reference genome with bowtie2 
# Note: Requires that reference is indexed (see README.md) 
# Note: Requires that input fastq files end in .fq.gz

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Global variables
TRIMMED_FOLDER="03_trimmed"
MAPPED_FOLDER="04_mapped"
REFERENCE="/home/ben/Documents/genomes/GCF_002872995.1_Otsh_v1.0_cds_from_genomic.fna"

# User variables
NUM_THREADS="4"

# Map reads and add RG
ls -1 $TRIMMED_FOLDER/*.fq.gz |
	sort -u |
	while read i
	do
	  echo $i
	  name=$(basename $i)
	  label=$(echo $name | perl -pe 's/\_trimmed\.fastq\.gz//')
	  ID="@RG\tID:${label}\tSM:${label}\tPL:Illumina"
	  bowtie2 --end-to-end -k 40 --threads $NUM_THREADS --rg-id $ID -x $REFERENCE -U $i -S $i.bowtie2.sam
	  samtools view -Sb $i.bowtie2.sam > $i.bowtie2.unsorted.bam
	  samtools sort -n -o $i.bowtie2.sorted.bam $i.bowtie2.unsorted.bam
      rm $i.bowtie2.sam
      rm $i.bowtie2.unsorted.bam
done

# clean up space
mv ./$TRIMMED_FOLDER/*.bam ./$MAPPED_FOLDER/

