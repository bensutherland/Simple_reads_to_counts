#!/bin/bash
# Use ref genome gtf to create known and unannotated output gtf for each sample
# Note: Requires sorted .bam 
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_log_files"

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Global variables
MAPPED_FOLDER="04_mapped"
REFERENCE_GFF="/home/ben/Documents/genomes/GCF_002872995.1_Otsh_v1.0_genomic.gff"
MERGELIST="00_archive/mergelist.txt"

# User variables
NUM_THREADS="3"

# Examine how the transcripts compare with the reference annotation
gffcompare -G -r $REFERENCE_GFF $MAPPED_FOLDER/stringtie_merged.gtf 

