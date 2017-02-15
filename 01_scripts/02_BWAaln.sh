#!/bin/bash
# Map trimmed reads to reference transcriptome with bwa mem
# Note: Requires that reference is indexed (see README.md) 

# Load modules
module load bwa/0.7.13
module load samtools/1.3

# Global variables
TRIMMED_FOLDER="03_trimmed"
MAPPED_FOLDER="04_mapped"
REFERENCE="/home/bensuth/00_resources/sfontinalis_contigs.fasta.gz"

# User variables
NUM_THREADS="5"


# map reads and add RG
ls -1 $TRIMMED_FOLDER/*.fastq.gz |
	sort -u |
	while read i
	do
	  echo $i
	  name=$(basename $i)
	  label=$(echo $name | perl -pe 's/.*(lib\d+[a-zA-Z]?).*/\1/')
	  ID="@RG\tID:${label}\tSM:${label}\tPL:Illumina"
	  bwa mem -t $NUM_THREADS -R ${ID} $REFERENCE $i > $i.sam
	  samtools view -Sb $i.sam > $i.unsorted.bam
	  samtools sort -o $i.sorted.bam $i.unsorted.bam
	  samtools index $i.sorted.bam
done

# clean up space
#rm ./$TRIMMED_FOLDER/*.sam ./$TRIMMED_FOLDER/*.unsorted.bam
mv ./$TRIMMED_FOLDER/*.bam ./$TRIMMED_FOLDER/*.sam ./$TRIMMED_FOLDER/*.bam.bai ./$MAPPED_FOLDER/
