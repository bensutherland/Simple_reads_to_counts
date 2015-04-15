#!/bin/bash
# use bwa to map trimmed reads to reference library

#set global variables
REFERENCE="sfontinalis_contigs.fasta"
READS="HI*_R1_trimmed.fastq.gz"
# OUTPUT="HI.2494.001.Index_12.lib06_R1_trimmed-out"

#Indexing reference library:  
gunzip -c ./$REFERENCE.gz > ./$REFERENCE

bwa index $REFERENCE
bwa mem -t 10 -R '@RG\tID:SAMPLENAME\tSM:SAMPLENAME\tPL:Illumina' $REFERENCE <samplenamehere> > samplename.sam

samtools import $REFERENCE.fai $OUTPUT.sam $OUTPUT.unsorted.bam
samtools sort $OUTPUT.unsorted.bam $OUTPUT
samtools index $OUTPUT.bam

# rm ./sfontinalis_contigs.fasta

