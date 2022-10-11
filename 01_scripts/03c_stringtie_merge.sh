#!/bin/bash
# Use ref genome GFF to create known and unannotated output GTF for each sample
# Note: Requires sorted .bam 

# Global variables
MAPPED_FOLDER="04_mapped"
REF_FOLDER="10_reference"
MERGELIST="00_archive/mergelist.txt"
NUM_THREADS="1" # seems to only use one thread

# User variables
REFERENCE_GFF="GCF_902806645.1_cgigas_uk_roslin_v1_genomic.gff"

# Run StringTie merge 
stringtie --merge -p $NUM_THREADS -G $REF_FOLDER/$REFERENCE_GFF \
	  -o $MAPPED_FOLDER/stringtie_merged.gtf \
	  $MERGELIST 

