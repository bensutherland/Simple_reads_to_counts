#!/bin/bash
# Map trimmed reads to reference transcriptome with bwa mem
# Note: Requires that reference is indexed (see README.md) 

# Global variables
TRIMMED_FOLDER="03_trimmed"
MAPPED_FOLDER="06_mapped"
REFERENCE="00_archive/reference_transcriptome/sfontinalis_contigs.fasta.gz"

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
	  samtools sort -o $i.sorted.bam $i.unsorted.bam
	  samtools index $i.sorted.bam
done

# clean up space
#rm ./$TRIMMED_FOLDER/*.sam ./$TRIMMED_FOLDER/*.unsorted.bam
mv ./$TRIMMED_FOLDER/*.bam ./$TRIMMED_FOLDER/*.bam.bai ./$MAPPED_FOLDER/
