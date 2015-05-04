#!/bin/bash
#first convert ordered BAM back into SAM

#set environment variables
MAPPED_FOLDER="07_mapped"
COUNT_FOLDER="09_GXlevels"

#convert sorted .bam back to .sam
ls -1 $MAPPED_FOLDER/*fastq.gz.bam | \
    sort -u | 
    while read i
    do
        echo $i
	samtools view -h "$i" > "$i".sam
    done

# mv $MAPPED_FOLDER/*_counts.txt $COUNT_FOLDER/
