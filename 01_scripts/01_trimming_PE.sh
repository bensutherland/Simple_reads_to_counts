#!/bin/bash
# Remove adapters and light quality trimming with Trimmomatic

# Global variables
RAW_FOLDER="02_raw_data"
TRIMMED_FOLDER="03_trimmed"
VECTORS="./00_archive/univec_trimmomatic.fasta"
TRIMMOMATIC_PROGRAM="/home/ben/Programs/trimmomatic-0.36.jar"

# User set variable
NUM_CPU="4"

# Filtering and trimming data with trimmomatic
ls -1 $RAW_FOLDER/*.fq.gz | 
    perl -pe 's/\_[12]\.fq\.gz//' |
    grep -vE "paired|single" |
    sort -u |
    while read i
    do
        echo "Trimming $i"
        java -Xmx35G -jar $TRIMMOMATIC_PROGRAM PE \
            -threads $NUM_CPU \
            -phred33 \
            "$i"_1.fq.gz \
            "$i"_2.fq.gz \
            "$i"_1.paired.fastq.gz \
            "$i"_1.single.fastq.gz \
            "$i"_2.paired.fastq.gz \
            "$i"_2.single.fastq.gz \
            ILLUMINACLIP:$VECTORS:2:30:10 \
            SLIDINGWINDOW:20:2 \
            LEADING:2 \
            TRAILING:2 \
            MINLEN:80
    done

# Move trimmed files to trimmed folder
mv $RAW_FOLDER/*single*.fastq.gz $TRIMMED_FOLDER
mv $RAW_FOLDER/*paired*.fastq.gz $TRIMMED_FOLDER
