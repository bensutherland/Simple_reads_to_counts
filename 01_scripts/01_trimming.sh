#!/bin/bash
# Remove adapters and light quality trimming with Trimmomatic
# Set location of trimmomatic.jar

# Global variables
RAW_FOLDER="02_raw_data"
TRIMMED_FOLDER="03_trimmed"
VECTORS="./00_archive/univec_trimmomatic.fasta"
TRIMMOMATIC_PROGRAM="/prg/trinityrnaseq/trinityrnaseq_r20140717/trinity-plugins/Trimmomatic-0.32/trimmomatic.jar"

# Filtering and trimming data with trimmomatic
ls -1 $RAW_FOLDER/*.fastq.gz | 
    sort -u |
    while read i
    do
        echo "Trimming $i"
        java -Xmx85G -jar $TRIMMOMATIC_PROGRAM SE \
            -threads 10 \
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
mv $RAW_FOLDER/*trimmed.fastq.gz $TRIMMED_FOLDER
