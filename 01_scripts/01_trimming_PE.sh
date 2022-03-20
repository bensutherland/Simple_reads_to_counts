#!/bin/bash
# Remove adapters and light quality trimming with Trimmomatic

# Global variables
RAW_FOLDER="02_raw_data"
TRIMMED_FOLDER="03_trimmed"
VECTORS="./00_archive/univec_trimmomatic.fasta"
TRIMMOMATIC_PROGRAM="/home/bsutherland/programs/Trimmomatic-0.39/trimmomatic-0.39.jar"

# User set variable
NUM_CPU="12"
EXTN="fastq.gz"

# Filtering and trimming data with trimmomatic (either fastq.gz or fq.gz)
ls -1 $RAW_FOLDER/{*.fq.gz,*.fastq.gz} | 

    # Remove extension
    perl -pe 's/[12]\.fq\.gz//' |
    perl -pe 's/[12]\.fastq\.gz//' |
    
    # Remove already-treated samples
    grep -vE "paired|single" |
    sort -u |
    while read i
    do
        echo "Trimming $i"
        java -Xmx35G -jar $TRIMMOMATIC_PROGRAM PE \
            -threads $NUM_CPU \
            -phred33 \
            "$i"1.$EXTN \
            "$i"2.$EXTN \
            "$i"1.paired.fq.gz \
            "$i"1.single.fq.gz \
            "$i"2.paired.fq.gz \
            "$i"2.single.fq.gz \
            ILLUMINACLIP:$VECTORS:2:30:10 \
            SLIDINGWINDOW:20:2 \
            LEADING:2 \
            TRAILING:2 \
            MINLEN:80
    done

# Move trimmed files to trimmed folder
mv $RAW_FOLDER/*single.fq.gz $TRIMMED_FOLDER
mv $RAW_FOLDER/*paired.fq.gz $TRIMMED_FOLDER
