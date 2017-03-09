#!/bin/bash
# Remove adapters and light quality trimming with Trimmomatic

# Global variables
RAW_FOLDER="02_raw_data"
TRIMMED_FOLDER="03_trimmed"
VECTORS="./00_archive/univec_trimmomatic.fasta"
TRIMMOMATIC_PROGRAM="/home/bensuth/programs/trimmomatic-0.36.jar"

# User set variable
NUM_CPU="5"


# Filtering and trimming data with trimmomatic
ls -1 $RAW_FOLDER/*.fastq.gz | 
    sort -u |
    while read i
    do
        echo "Trimming $i"
        java -Xmx75G -jar $TRIMMOMATIC_PROGRAM SE \
            -threads $NUM_CPU \
            -phred33 \
            "$i" \
            "${i%.fastq.gz}"_trimmed.fastq.gz \
            ILLUMINACLIP:$VECTORS:2:30:10 \
            SLIDINGWINDOW:20:2 \
            LEADING:2 \
            TRAILING:2 \
            MINLEN:80
    done

# Move trimmed files to trimmed folder
mv $RAW_FOLDER/*_trimmed.fastq.gz $TRIMMED_FOLDER
