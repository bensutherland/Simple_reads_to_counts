#!/bin/bash
# use bwa mem to map trimmed reads to reference library

###NOTE THIS SCRIPT REQUIRES THAT YOUR REFERENCE IS ALREADY INDEXED (bwa index REFERENCE)

# global variables (note: point to REFERENCE)
TRIMMED_FOLDER="03_trimmed"
MAPPED_FOLDER="06_mapped"
REFERENCE="05_trinity_output/sfontinalis_contigs.fasta"

# map reads and add RG
ls -1 $TRIMMED_FOLDER/*.fastq.gz |
	sort -u |
	while read i
	do
	  echo $i
	  name=$(basename $i)
	  label=$(echo $name | perl -pe 's/.*(lib\d+[a-zA-Z]?).*/\1/')
	  ID="@RG\tID:${label}\tSM:${label}\tPL:Illumina"
	  bwa mem -t 10 -R ${ID} $REFERENCE $i > $i.sam
	  samtools view -Sb $i.sam > $i.unsorted.bam
	  samtools sort $i.unsorted.bam $i
	  samtools index $i.bam
done

# clean up space
rm ./$TRIMMED_FOLDER/*.sam ./$TRIMMED_FOLDER/*.unsorted.bam
mv ./$TRIMMED_FOLDER/*.bam ./$TRIMMED_FOLDER/*.bam.bai ./$MAPPED_FOLDER/
