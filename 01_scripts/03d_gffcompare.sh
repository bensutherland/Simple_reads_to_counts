#!/bin/bash
# Generate statistics based on your de novo and reference transcripts
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Global variables
MAPPED_FOLDER="04_mapped"
REFERENCE_GFF="/scratch2/bsutherland/ref_genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.gff"
MERGELIST="00_archive/mergelist.txt"

# User variables
NUM_THREADS="3"

# Examine how the transcripts compare with the reference annotation
gffcompare -G -r $REFERENCE_GFF $MAPPED_FOLDER/stringtie_merged.gtf 

