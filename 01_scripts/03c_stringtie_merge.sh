#!/bin/bash
# Use ref genome gtf to create known and unannotated output gtf for each sample
# Note: Requires sorted .bam 
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Variables
MAPPED_FOLDER="04_mapped"
REFERENCE_GFF="/scratch2/bsutherland/ref_genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.gff"
MERGELIST="00_archive/mergelist.txt"
NUM_THREADS="1" # seems to only use one thread

# Run Stringtie merge 
stringtie --merge -p $NUM_THREADS -G $REFERENCE_GFF -o $MAPPED_FOLDER/stringtie_merged.gtf $MERGELIST 

