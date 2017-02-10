#!/bin/bash
# Remove adapters and light quality trimming with Trimmomatic

# Load module 
module load trimmomatic/0.33 

# Global variables
RAW_FOLDER="02_raw_data"
TRIMMED_FOLDER="03_trimmed"
VECTORS="./00_archive/univec_trimmomatic.fasta"
TRIMMOMATIC_PROGRAM="trimmomatic-0.33.jar"

# User set variable
NUM_CPU="10"


# Filtering and trimming data with trimmomatic
ls -1 $RAW_FOLDER/*.fastq.gz | 
    sort -u |
    while read i
    do
        echo "Trimming $i"
        $TRIMMOMATIC_PROGRAM SE \
            -threads $NUM_CPU \
            -phred33 \
            "$i" \
            "${i%R1.fastq.gz}"R1_trimmed.fastq.gz \
            ILLUMINACLIP:$VECTORS:2:30:10 \
            SLIDINGWINDOW:20:2 \
            LEADING:2 \
            TRAILING:2 \
            MINLEN:80
    done

# Move trimmed files to trimmed folder
mv $RAW_FOLDER/*R1_trimmed.fastq.gz $TRIMMED_FOLDER
