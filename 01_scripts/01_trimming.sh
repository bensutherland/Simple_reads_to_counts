#!/bin/bash
# Cleaning and trimming fastq.gz read files with Trimmomatic
# point to the trimmomatic.jar in Global variables

# Global variables
RAW_FOLDER="02_raw_data"
TRIMMED_FOLDER="03_trimmed"
VECTORS="./00_archive/univec_trimmomatic.fasta"
TRIMMOMATIC_PROGRAM="/prg/trinityrnaseq/trinityrnaseq_r20140717/trinity-plugins/Trimmomatic-0.32/trimmomatic.jar"

# Filtering and trimming data with trimmomatic
ls -1 $RAW_FOLDER/*.fastq.gz | \
    sort -u | \ 
    while read i
    do
        echo $i
        java -Xmx70G -jar $TRIMMOMATIC_PROGRAM SE \
            -threads 10 \
            -phred33 \
            "$i" \
            "${i%R1.fastq.gz}"R1_trimmed.fastq.gz \
            ILLUMINACLIP:$VECTORS:2:30:10 \
            SLIDINGWINDOW:20:30 \
            LEADING:20 \
            TRAILING:20 \
            MINLEN:80
    done

# Transfer trimmed files to $TRIMMED_FOLDER folder
mv $RAW_FOLDER/*trimmed* $TRIMMED_FOLDER
